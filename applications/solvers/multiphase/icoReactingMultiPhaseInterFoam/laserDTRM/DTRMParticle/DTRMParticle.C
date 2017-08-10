/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenCFD Ltd
     \\/     M anipulation  |
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

#include "DTRMParticle.H"
#include "constants.H"
#include "physicoChemicalConstants.H"

using namespace Foam::constant;
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DTRMParticle::DTRMParticle
(
    const polyMesh& mesh,
    const vector& position,
    const vector& targetPosition,
    const scalar I,
    const label cellI,
    const scalar dA,
    const label reflectedId,
    const scalar Imin,
    bool doCellFacePt
)
:
    particle(mesh, position, cellI, doCellFacePt),
    p0_(position),
    p1_(targetPosition),
    I0_(I),
    I_(I),
    dA_(dA),
    reflectedId_(reflectedId),
    Imin_(Imin)
{}


Foam::DTRMParticle::DTRMParticle(const DTRMParticle& p)
:
    particle(p),
    p0_(p.p0_),
    p1_(p.p1_),
    I0_(p.I0_),
    I_(p.I_),
    dA_(p.dA_),
    reflectedId_(p.reflectedId_),
    Imin_(p.Imin_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::DTRMParticle::move
(
    trackingData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    const polyBoundaryMesh& pbMesh = mesh_.boundaryMesh();

    while (td.keepParticle && !td.switchProcessor)
    {
        point p0 = position();
        label cell0 = cell();

        scalar dt = trackToFace(p1_, td);

        // Consider the cell between f0(start of tracking) and f1
        label celli = cell();

        const vector dsv = position() - p0;
        const scalar ds = mag(dsv);

        // Boltzman constant
        const scalar sigma = physicoChemical::sigma.value();
        if
        (
            (!td.relfectedCells()[celli] > 0 && reflectedId_ == 0)
         || reflectedId_ > 0
        )
        {
            scalar a = td.aInterp().interpolate(position(), cell0);
            scalar e = td.eInterp().interpolate(position(), cell0);
            scalar E = td.EInterp().interpolate(position(), cell0);
            scalar T = td.TInterp().interpolate(position(), cell0);

            const scalar I1 =
            (
                I_
                + ds*(e*sigma*pow4(T)/mathematical::pi + E)
            ) / (1 + ds*a);

            td.Q(cell0) += (I_ - I1)*dA_;

            I_ = I1;

            if ((I_ <= 0.01*I0_) || (I_ < Imin_))
            {
                break;
            }
        }
        else
        {
            scalar rho(0);
            // Create a new reflected particle when the particles is not
            // transmissive and larger than an absolute I
            if (reflectedId_ == 0 && I_ > Imin_)
            {
                vector pDir = dsv/ds;

                cellPointWeight cpw(mesh_, position(), celli, face());
                //vector nHat = td.nHatCells()[celli];
                vector nHat = td.nHatInterp().interpolate(cpw);

                nHat /= mag(nHat);
                scalar cosTheta(-pDir & nHat);
                // Only new incoming rays
                if (cosTheta > SMALL)
                {
                    vector newDir = td.reflection().R(pDir, nHat);

                    //scalar theta = acos(-pDir & nHat);

                    // reflectivity
                    rho = min(max(td.reflection().rho(cosTheta), 0.0), 0.98);

                    scalar delaM = sqrt(mesh_.cellVolumes()[cell0]);

                    DTRMParticle* pPtr = new DTRMParticle
                    (
                        mesh_,
                        position() - pDir*0.1*delaM,
                        position() + newDir*mesh_.bounds().mag(),
                        I_*rho,
                        cell0,
                        dA_,
                        reflectedId_,
                        Imin_,
                        true
                    );
                    // Add to cloud
                    td.cloud().addParticle(pPtr);
                }
            }

            reflectedId_++;

            const point p0 = position();

            // Try to locate this particle across the reflecting surface in
            // a pure phase face
            scalar dt = trackToFace(p1_, td);
            const scalar ds = mag(position() - p0);

            scalar a = td.aInterp().interpolate(position(), celli);
            scalar e = td.eInterp().interpolate(position(), celli);
            scalar E = td.EInterp().interpolate(position(), celli);
            scalar T = td.TInterp().interpolate(position(), celli);

            // Left intensity after reflection
            const scalar Itran = I_*(1.0 - rho);
            const scalar I1 =
            (
                Itran
              + ds*(e*sigma*pow4(T)/mathematical::pi + E)
            ) / (1 + ds*a);

            td.Q(celli) += (Itran - I1)*dA_;

            I_ = I1;

            if (I_ <= 0.01*I0_ || I_ < Imin_)
            {
                break;
            }
        }

        if (onBoundary() && td.keepParticle)
        {
            if (isA<processorPolyPatch>(pbMesh[patch(face())]))
            {
                td.switchProcessor = true;
            }
        }
    }

    return td.keepParticle;
}


bool Foam::DTRMParticle::hitPatch
(
    const polyPatch&,
    particle::TrackingData<Cloud<DTRMParticle>>&,
    const label,
    const scalar,
    const tetIndices&
)
{
    return false;
}


void Foam::DTRMParticle::hitProcessorPatch
(
    const processorPolyPatch&,
    particle::TrackingData<Cloud<DTRMParticle>>& td
)
{
    td.switchProcessor = true;
}


void Foam::DTRMParticle::hitWallPatch
(
    const wallPolyPatch& wpp,
    particle::TrackingData<Cloud<DTRMParticle>>& td,
    const tetIndices& tetIs
)
{
    td.keepParticle = false;
}


void Foam::DTRMParticle::hitPatch
(
    const polyPatch&,
    particle::TrackingData<Cloud<DTRMParticle>>& td
)
{
    td.keepParticle = false;
}


// ************************************************************************* //
