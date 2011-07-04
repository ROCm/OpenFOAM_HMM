/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "monoatomic.H"
#include "Random.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(monoatomic, 0);
    defineTemplateTypeNameAndDebug(Cloud<monoatomic>, 0);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::monoatomic::move
(
    monoatomic::trackingData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    if (special_ != SPECIAL_FROZEN)
    {
        return td.keepParticle;
    }

    const constantProperties& constProps(td.cloud().constProps(id_));

    if (td.part() == trackingData::tpFirstVelocityHalfStep)
    {
        // First leapfrog velocity adjust part, required before tracking+force
        // part

        v_ += 0.5*trackTime*a_;
    }
    else if (td.part() == trackingData::tpLinearTrack)
    {
        // Leapfrog tracking part

        scalar tEnd = (1.0 - stepFraction())*trackTime;
        scalar dtMax = tEnd;

        while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
        {
            // set the lagrangian time-step
            scalar dt = min(dtMax, tEnd);

            dt *= trackToFace(position() + dt*v_, td);

            tEnd -= dt;
            stepFraction() = 1.0 - tEnd/trackTime;
        }

        setSitePositions(constProps);
    }
    else if (td.part() == trackingData::tpSecondVelocityHalfStep)
    {
        // Second leapfrog velocity adjust part, required after tracking+force
        // part

        a_ = siteForces_[0]/constProps.mass();

        v_ += 0.5*trackTime*a_;
    }
    else if (td.part() != trackingData::tpRotationalTrack)
    {
        FatalErrorIn("monoatomic::move(trackingData&, const scalar)") << nl
            << td.part() << " is an invalid part of the integration method."
            << abort(FatalError);
    }

    return td.keepParticle;
}


void Foam::monoatomic::transformProperties(const tensor& T)
{
    particle::transformProperties(T);

    v_ = transform(T, v_);

    a_ = transform(T, a_);

    rf_ = transform(T, rf_);

    sitePositions_[0] = position_ + (T & (sitePositions_[0] - position_));

    siteForces_[0] = T & siteForces_[0];
}


void Foam::monoatomic::transformProperties(const vector& separation)
{
    particle::transformProperties(separation);

    if (special_ == SPECIAL_TETHERED)
    {
        specialPosition_ += separation;
    }

    sitePositions_[0] += separation;
}


void Foam::monoatomic::setSitePositions(const constantProperties& constProps)
{
    sitePositions_[0] = position_;
}


void Foam::monoatomic::setSiteSizes(label size)
{
    // Nothing required, size controlled internally
}


bool Foam::monoatomic::hitPatch
(
    const polyPatch&,
    trackingData&,
    const label,
    const scalar,
    const tetIndices&
)
{
    return false;
}


void Foam::monoatomic::hitProcessorPatch
(
    const processorPolyPatch&,
    trackingData& td
)
{
    td.switchProcessor = true;
}


void Foam::monoatomic::hitWallPatch
(
    const wallPolyPatch& wpp,
    trackingData& td,
    const tetIndices& tetIs
)
{
    // Use of the normal from tetIs is not required as
    // hasWallImpactDistance is false.
    vector nw = normal();
    nw /= mag(nw);

    scalar vn = v_ & nw;

    // Specular reflection
    if (vn > 0)
    {
        v_ -= 2*vn*nw;
    }
}


void Foam::monoatomic::hitPatch
(
    const polyPatch&,
    trackingData& td
)
{
    td.keepParticle = false;
}


// ************************************************************************* //
