/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "FreeStream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::FreeStream<CloudType>::FreeStream
(
    const dictionary& dict,
    CloudType& cloud
)
:
    InflowBoundaryModel<CloudType>(dict, cloud, typeName),
    patchIndex_(),
    temperature_(readScalar(this->coeffDict().lookup("temperature"))),
    velocity_(this->coeffDict().lookup("velocity")),
    moleculeTypeIds_(),
    numberDensities_(),
    particleFluxAccumulators_()
{
    word patchName = this->coeffDict().lookup("patch");

    patchIndex_ = cloud.mesh().boundaryMesh().findPatchID(patchName);

    const polyPatch& patch = cloud.mesh().boundaryMesh()[patchIndex_];

    if (patchIndex_ == -1)
    {
        FatalErrorIn("Foam::DsmcCloud<ParcelType>::initialise")
        << "patch " << patchName << " not found." << nl
            << abort(FatalError);
    }

    const dictionary& numberDensitiesDict
    (
        this->coeffDict().subDict("numberDensities")
    );

    List<word> molecules(numberDensitiesDict.toc());

    numberDensities_.setSize(molecules.size());

    moleculeTypeIds_.setSize(molecules.size());

    forAll(molecules, i)
    {
        numberDensities_[i] = readScalar
        (
            numberDensitiesDict.lookup(molecules[i])
        );

        moleculeTypeIds_[i] = findIndex(cloud.typeIdList(), molecules[i]);

        if (moleculeTypeIds_[i] == -1)
        {
            FatalErrorIn("Foam::DsmcCloud<ParcelType>::initialise")
                << "typeId " << molecules[i] << "not defined in cloud." << nl
                << abort(FatalError);
        }
    }

    numberDensities_ /= cloud.nParticle();

    particleFluxAccumulators_.setSize
    (
        molecules.size(),
        Field<scalar>(patch.size(), 0)
    );
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::FreeStream<CloudType>::~FreeStream()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class CloudType>
void Foam::FreeStream<CloudType>::inflow()
{
    CloudType& cloud(this->owner());

    const polyMesh& mesh(cloud.mesh());

    const scalar deltaT = mesh.time().deltaT().value();

    Random& rndGen(cloud.rndGen());

    const polyPatch& patch = mesh.boundaryMesh()[patchIndex_];

    label particlesInserted = 0;

    // Add mass to the accumulators.  negative face area dotted with the
    // velocity to point flux into the domain.

    forAll(particleFluxAccumulators_, i)
    {
        particleFluxAccumulators_[i] +=
           -patch.faceAreas() & (velocity_*numberDensities_[i]*deltaT);
    }

    forAll(patch, f)
    {
        // Loop over all faces as the outer loop to avoid calculating
        // geometrical properties multiple times.

        labelList faceVertices = patch[f];

        label nVertices = faceVertices.size();

        label globalFaceIndex = f + patch.start();

        label cell = mesh.faceOwner()[globalFaceIndex];

        const vector& fC = patch.faceCentres()[f];

        scalar fA = mag(patch.faceAreas()[f]);

        // Cummulative triangle area fractions
        List<scalar> cTriAFracs(nVertices);

        for (label v = 0; v < nVertices - 1; v++)
        {
            const point& vA = mesh.points()[faceVertices[v]];

            const point& vB = mesh.points()[faceVertices[(v + 1)]];

            cTriAFracs[v] =
                0.5*mag((vA - fC)^(vB - fC))/fA
              + cTriAFracs[max((v - 1), 0)];
        }

        // Force the last area fraction value to 1.0 to avoid any
        // rounding/non-flat face errors giving a value < 1.0
        cTriAFracs[nVertices - 1] = 1.0;

        // Normal unit vector *negative* so normal is pointing into the
        // domain
        vector nw = patch.faceAreas()[f];
        nw /= -mag(nw);

        // Wall tangential unit vector. Use the direction between the
        // face centre and the first vertex in the list
        vector tw1 = fC - (mesh.points()[faceVertices[0]]);
        tw1 /= mag(tw1);

        // Other tangential unit vector.  Rescaling in case face is not
        // flat and nw and tw1 aren't perfectly orthogonal
        vector tw2 = nw ^ tw1;
        tw2 /= mag(tw2);

        forAll(particleFluxAccumulators_, i)
        {
            scalar& faceAccumulator = particleFluxAccumulators_[i][f];

            // Number of particles to insert
            label nI = max(label(faceAccumulator), 0);

            faceAccumulator -= nI;

            label typeId = moleculeTypeIds_[i];

            scalar mass = cloud.constProps(typeId).mass();

            for (label n = 0; n < nI; n++)
            {
                // Choose a triangle to insert on, based on their relative area

                scalar triSelection = rndGen.scalar01();

                // Selected triangle
                label sTri = -1;

                forAll(cTriAFracs, tri)
                {
                    sTri = tri;

                    if (cTriAFracs[tri] >= triSelection)
                    {
                        break;
                    }
                }

                // Randomly distribute the points on the triangle, using the
                // method from:
                // Generating Random Points in Triangles
                // by Greg Turk
                // from "Graphics Gems", Academic Press, 1990
                // http://tog.acm.org/GraphicsGems/gems/TriPoints.c

                const point& A = fC;
                const point& B = mesh.points()[faceVertices[sTri]];
                const point& C =
                    mesh.points()[faceVertices[(sTri + 1) % nVertices]];

                scalar s = rndGen.scalar01();
                scalar t = sqrt(rndGen.scalar01());

                point p = (1 - t)*A + (1 - s)*t*B + s*t*C;

                vector U =
                    sqrt(CloudType::kb*temperature_/mass)
                   *(
                        rndGen.GaussNormal()*tw1
                      + rndGen.GaussNormal()*tw2
                      - sqrt(-2.0*log(max(1 - rndGen.scalar01(),VSMALL)))*nw
                    );

                U += velocity_;

                scalar Ei = cloud.equipartitionInternalEnergy
                (
                    temperature_,
                    cloud.constProps(typeId).internalDegreesOfFreedom()
                );

                cloud.addNewParcel
                (
                    p,
                    U,
                    Ei,
                    cell,
                    typeId
                );

                particlesInserted++;
            }
        }
    }

    reduce(particlesInserted, sumOp<label>());

    Info<< "    Particles inserted              = "
        << particlesInserted << endl;

}


// ************************************************************************* //
