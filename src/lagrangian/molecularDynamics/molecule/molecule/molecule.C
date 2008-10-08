/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "moleculeCloud.H"
#include "molecule.H"
#include "Random.H"
#include "Time.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tensor Foam::molecule::rotationTensor(scalar deltaT) const
{
    scalar phi1 = 0.5*deltaT*omega_.x();

    tensor U1
    (
        tensor
        (
            1, 0, 0,
            0, Foam::cos(phi1), -Foam::sin(phi1),
            0, Foam::sin(phi1), Foam::cos(phi1)
        )
    );

    scalar phi2 = 0.5*deltaT*omega_.y();

    tensor U2
    (
        tensor
        (
            Foam::cos(phi2), 0, Foam::sin(phi2),
            0, 1, 0,
            -Foam::sin(phi2), 0, Foam::cos(phi2)
        )
    );

    scalar phi3 = deltaT*omega_.z();

    tensor U3
    (
        tensor
        (
            Foam::cos(phi3), -Foam::sin(phi3), 0,
            Foam::sin(phi3), Foam::cos(phi3), 0,
            0, 0, 1
        )
    );

    return (U1 & U2 & U3 & U2 & U1);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::molecule::trackData::trackData
(
    moleculeCloud& molCloud,
    label part
)
:
    Particle<molecule>::trackData(molCloud),
    molCloud_(molCloud),
    part_(part)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::molecule::move(molecule::trackData& td)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    scalar deltaT = cloud().pMesh().time().deltaT().value();
    scalar tEnd = (1.0 - stepFraction())*deltaT;
    scalar dtMax = tEnd;

    if (td.part() == 0)
    {
        // First leapfrog velocity adjust part, required before tracking+force
        // part

        v_ += 0.5*deltaT*a_;

        omega_ += 0.5*deltaT*alpha_;
    }
    else if (td.part() == 1)
    {
        // Leapfrog tracking part

        while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
        {
            // set the lagrangian time-step
            scalar dt = min(dtMax, tEnd);

            dt *= trackToFace(position() + dt*v_, td);

            tEnd -= dt;
            stepFraction() = 1.0 - tEnd/deltaT;
        }
    }
    else if (td.part() == 2)
    {
        // Leapfrog orientation adjustment, carried out before force calculation
        // but after tracking stage, i.e. rotation carried once linear motion
        // complete.

        R_ = R_ & rotationTensor(deltaT);

        setSitePositions(td.molCloud().constProps(id_));
    }
    else if (td.part() == 3)
    {
        // Second leapfrog velocity adjust part, required after tracking+force
        // part

        const diagTensor& I(td.molCloud().constProps(id_).momentOfInertia());

        scalar m = td.molCloud().constProps(id_).mass();

        a_ = vector::zero;

        vector tau(vector::zero);

        forAll(siteForces_, s)
        {
            const vector& f = siteForces_[s];

            a_ += f/m;

            tau += (sitePositions_[s] - position_) ^ f;
        }

        alpha_ = R_ & inv(I) & R_.T() & tau;

        v_ += 0.5*deltaT*a_;

        omega_ += 0.5*deltaT*alpha_;
    }
    else
    {
        FatalErrorIn("molecule::move(molecule::trackData& td)") << nl
            << td.part()
            << " is an invalid part of the integration method."
            << abort(FatalError);
    }

    return td.keepParticle;
}


void Foam::molecule::transformProperties(const tensor& T)
{}


void Foam::molecule::transformProperties(const vector& separation)
{
    if (special_ == SPECIAL_TETHERED)
    {
        specialPosition_ += separation;
    }
}


void Foam::molecule::setSitePositions(const constantProperties& constProps)
{
    sitePositions_ = position_ + (R_ & constProps.siteReferencePositions());
}


void Foam::molecule::hitProcessorPatch
(
    const processorPolyPatch&,
    molecule::trackData& td
)
{
    td.switchProcessor = true;
}


void Foam::molecule::hitProcessorPatch
(
    const processorPolyPatch&,
    int&
)
{}


void Foam::molecule::hitWallPatch
(
    const wallPolyPatch& wpp,
    molecule::trackData& td
)
{
    vector nw = wpp.faceAreas()[wpp.whichFace(face())];
    nw /= mag(nw);

    scalar vn = v_ & nw;
//     vector vt = v_ - vn*nw;

//     Random rand(clock::getTime());

//     scalar tmac = 0.8;

//     scalar wallTemp = 2.5;

//     if (rand.scalar01() < tmac)
//     {
//         // Diffuse reflection
//
//         vector tw1 = vt/mag(vt);
//
//         vector tw2 = nw ^ tw1;
//
//         V_ = sqrt(wallTemp/mass_)*rand.GaussNormal()*tw1
//                 + sqrt(wallTemp/mass_)*rand.GaussNormal()*tw2
//                 - mag(sqrt(wallTemp/mass_)*rand.GaussNormal())*nw;
//     }

//     else
//     {
        // Specular reflection

        if (vn > 0)
        {
            v_ -= 2*vn*nw;
        }

//     }

}


void Foam::molecule::hitWallPatch
(
    const wallPolyPatch&,
    int&
)
{}


void Foam::molecule::hitPatch
(
    const polyPatch&,
    molecule::trackData& td
)
{
    td.keepParticle = false;
}


void Foam::molecule::hitPatch
(
    const polyPatch&,
    int&
)
{}


// ************************************************************************* //
