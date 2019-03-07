/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "cyclicACMIFvPatch.H"
#include "fvMesh.H"
#include "transform.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicACMIFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, cyclicACMIFvPatch, polyPatch);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::cyclicACMIFvPatch::resetPatchAreas(const fvPatch& fvp) const
{
    const_cast<vectorField&>(fvp.Sf()) = fvp.patch().faceAreas();
    const_cast<scalarField&>(fvp.magSf()) = mag(fvp.patch().faceAreas());

    DebugPout
        << fvp.patch().name() << " area:" << sum(fvp.magSf()) << endl;
}


void Foam::cyclicACMIFvPatch::makeWeights(scalarField& w) const
{
    if (coupled())
    {
        const cyclicACMIFvPatch& nbrPatch = neighbFvPatch();
        const scalarField deltas(nf() & coupledFvPatch::delta());

        // These deltas are of the cyclic part alone - they are
        // not affected by the amount of overlap with the nonOverlapPatch
        scalarField nbrDeltas
        (
            interpolate
            (
                nbrPatch.nf() & nbrPatch.coupledFvPatch::delta()
            )
        );

        const scalar tol = cyclicACMIPolyPatch::tolerance();


        forAll(deltas, facei)
        {
            scalar di = mag(deltas[facei]);
            scalar dni = mag(nbrDeltas[facei]);

            if (dni < tol)
            {
                // Avoid zero weights on disconnected faces. This value
                // will be weighted with the (zero) face area so will not
                // influence calculations.
                w[facei] = 1.0;
            }
            else
            {
                w[facei] = dni/(di + dni);
            }
        }
    }
    else
    {
        // Behave as uncoupled patch
        fvPatch::makeWeights(w);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cyclicACMIFvPatch::coupled() const
{
    return Pstream::parRun() || (this->size() && neighbFvPatch().size());
}


Foam::tmp<Foam::vectorField> Foam::cyclicACMIFvPatch::delta() const
{
    if (coupled())
    {
        const cyclicACMIFvPatch& nbrPatch = neighbFvPatch();

        const vectorField patchD(coupledFvPatch::delta());

        vectorField nbrPatchD(interpolate(nbrPatch.coupledFvPatch::delta()));

        tmp<vectorField> tpdv(new vectorField(patchD.size()));
        vectorField& pdv = tpdv.ref();

        // do the transformation if necessary
        if (parallel())
        {
            forAll(patchD, facei)
            {
                const vector& ddi = patchD[facei];
                const vector& dni = nbrPatchD[facei];

                pdv[facei] = ddi - dni;
            }
        }
        else
        {
            forAll(patchD, facei)
            {
                const vector& ddi = patchD[facei];
                const vector& dni = nbrPatchD[facei];

                pdv[facei] = ddi - transform(forwardT()[0], dni);
            }
        }

        return tpdv;
    }
    else
    {
        return coupledFvPatch::delta();
    }
}


Foam::tmp<Foam::labelField> Foam::cyclicACMIFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}


Foam::tmp<Foam::labelField> Foam::cyclicACMIFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    return neighbFvPatch().patchInternalField(iF);
}


void Foam::cyclicACMIFvPatch::movePoints()
{
    if (!cyclicACMIPolyPatch_.owner())
    {
        return;
    }

    // Set the patch face areas to be consistent with the changes made at the
    // polyPatch level

    const fvPatch& nonOverlapPatch = this->nonOverlapPatch();
    const cyclicACMIFvPatch& nbrACMI = neighbPatch();
    const fvPatch& nbrNonOverlapPatch = nbrACMI.nonOverlapPatch();

    resetPatchAreas(*this);
    resetPatchAreas(nonOverlapPatch);
    resetPatchAreas(nbrACMI);
    resetPatchAreas(nbrNonOverlapPatch);

    // Scale the mesh flux

    const labelListList& newSrcAddr = AMI().srcAddress();
    const labelListList& newTgtAddr = AMI().tgtAddress();

    const fvMesh& mesh = boundaryMesh().mesh();
    surfaceScalarField& meshPhi = const_cast<fvMesh&>(mesh).setPhi();
    surfaceScalarField::Boundary& meshPhiBf = meshPhi.boundaryFieldRef();

    // Note: phip and phiNonOverlapp will be different sizes if new faces
    // have been added
    scalarField& phip = meshPhiBf[cyclicACMIPolyPatch_.index()];
    scalarField& phiNonOverlapp =
        meshPhiBf[nonOverlapPatch.patch().index()];

    forAll(phip, facei)
    {
        if (newSrcAddr[facei].empty())
        {
            // AMI patch with no connection to other coupled faces
            phip[facei] = 0.0;
        }
        else
        {
            // Scale the mesh flux according to the area fraction
            const face& fAMI = cyclicACMIPolyPatch_.localFaces()[facei];

            // Note: using raw point locations to calculate the geometric
            // area - faces areas are currently scaled (decoupled from
            // mesh points)
            const scalar geomArea = fAMI.mag(cyclicACMIPolyPatch_.localPoints());
            phip[facei] *= magSf()[facei]/geomArea;
        }
    }

    forAll(phiNonOverlapp, facei)
    {
        const scalar w = 1.0 - cyclicACMIPolyPatch_.srcMask()[facei];
        phiNonOverlapp[facei] *= w;
    }

    scalarField& nbrPhip = meshPhiBf[nbrACMI.patch().index()];
    scalarField& nbrPhiNonOverlapp =
        meshPhiBf[nbrNonOverlapPatch.patch().index()];

    forAll(nbrPhip, facei)
    {
        if (newTgtAddr[facei].empty())
        {
            nbrPhip[facei] = 0.0;
        }
        else
        {
            const face& fAMI = nbrACMI.patch().localFaces()[facei];

            // Note: using raw point locations to calculate the geometric
            // area - faces areas are currently scaled (decoupled from
            // mesh points)
            const scalar geomArea = fAMI.mag(nbrACMI.patch().localPoints());
            nbrPhip[facei] *= nbrACMI.magSf()[facei]/geomArea;
        }
    }

    forAll(nbrPhiNonOverlapp, facei)
    {
        const scalar w = 1.0 - cyclicACMIPolyPatch_.tgtMask()[facei];
        nbrPhiNonOverlapp[facei] *= w;
    }
}


// ************************************************************************* //
