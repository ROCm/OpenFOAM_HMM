/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "wallBoundedStreamLineParticle.H"
#include "vectorFieldIOField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::vector Foam::wallBoundedStreamLineParticle::interpolateFields
(
    const trackingData& td,
    const point& position,
    const label celli,
    const label facei
)
{
    if (celli == -1)
    {
        FatalErrorInFunction
            << "Cell:" << celli << abort(FatalError);
    }

    const vector U =
        td.vvInterp_[td.UIndex_].interpolate(position, celli, facei);

    // Check if at different position
    if
    (
       !sampledPositions_.size()
     || magSqr(sampledPositions_.last() - position) > Foam::sqr(SMALL)
    )
    {
        // Store the location
        sampledPositions_.append(position);

        // Store the scalar fields
        sampledScalars_.setSize(td.vsInterp_.size());
        forAll(td.vsInterp_, scalari)
        {
            sampledScalars_[scalari].append
            (
                td.vsInterp_[scalari].interpolate(position, celli, facei)
            );
        }

        // Store the vector fields
        sampledVectors_.setSize(td.vvInterp_.size());
        forAll(td.vvInterp_, vectori)
        {
            vector positionU;
            if (vectori == td.UIndex_)
            {
                positionU = U;
            }
            else
            {
                positionU =
                    td.vvInterp_[vectori].interpolate(position, celli, facei);
            }
            sampledVectors_[vectori].append(positionU);
        }
    }

    return U;
}


Foam::vector Foam::wallBoundedStreamLineParticle::sample
(
    trackingData& td
)
{
    vector U = interpolateFields(td, localPosition_, cell(), face());

    if (!td.trackForward_)
    {
        U = -U;
    }

    scalar magU = mag(U);

    if (magU < SMALL)
    {
        // Stagnant particle. Might as well stop
        lifeTime_ = 0;
        return vector::zero;
    }
    else
    {
        return U/magU;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoundedStreamLineParticle::wallBoundedStreamLineParticle
(
    const polyMesh& mesh,
    const point& position,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const label meshEdgeStart,
    const label diagEdge,
    const label lifeTime
)
:
    wallBoundedParticle
    (
        mesh,
        position,
        celli,
        tetFacei,
        tetPti,
        meshEdgeStart,
        diagEdge
    ),
    lifeTime_(lifeTime)
{}


Foam::wallBoundedStreamLineParticle::wallBoundedStreamLineParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    wallBoundedParticle(mesh, is, readFields, newFormat)
{
    if (readFields)
    {
        List<scalarList> sampledScalars;
        List<vectorList> sampledVectors;

        is  >> lifeTime_
            >> sampledPositions_ >> sampledScalars >> sampledVectors;

        sampledScalars_.setSize(sampledScalars.size());
        forAll(sampledScalars, i)
        {
            sampledScalars_[i].transfer(sampledScalars[i]);
        }
        sampledVectors_.setSize(sampledVectors.size());
        forAll(sampledVectors, i)
        {
            sampledVectors_[i].transfer(sampledVectors[i]);
        }
    }

    is.check(FUNCTION_NAME);
}


Foam::wallBoundedStreamLineParticle::wallBoundedStreamLineParticle
(
    const wallBoundedStreamLineParticle& p
)
:
    wallBoundedParticle(p),
    lifeTime_(p.lifeTime_),
    sampledPositions_(p.sampledPositions_),
    sampledScalars_(p.sampledScalars_),
    sampledVectors_(p.sampledVectors_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wallBoundedStreamLineParticle::readFields
(
    Cloud<wallBoundedStreamLineParticle>& c
)
{
    if (!c.size())
    {
        return;
    }

    wallBoundedParticle::readFields(c);

    IOField<label> lifeTime
    (
        c.fieldIOobject("lifeTime", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, lifeTime);

    vectorFieldIOField sampledPositions
    (
        c.fieldIOobject("sampledPositions", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, sampledPositions);

    label i = 0;
    forAllIter(Cloud<wallBoundedStreamLineParticle>, c, iter)
    {
        iter().lifeTime_ = lifeTime[i];
        iter().sampledPositions_.transfer(sampledPositions[i]);
        i++;
    }
}


void Foam::wallBoundedStreamLineParticle::writeFields
(
    const Cloud<wallBoundedStreamLineParticle>& c
)
{
    wallBoundedParticle::writeFields(c);

    label np =  c.size();

    IOField<label> lifeTime
    (
        c.fieldIOobject("lifeTime", IOobject::NO_READ),
        np
    );
    vectorFieldIOField sampledPositions
    (
        c.fieldIOobject("sampledPositions", IOobject::NO_READ),
        np
    );

    label i = 0;
    forAllConstIter(Cloud<wallBoundedStreamLineParticle>, c, iter)
    {
        lifeTime[i] = iter().lifeTime_;
        sampledPositions[i] = iter().sampledPositions_;
        i++;
    }

    lifeTime.write();
    sampledPositions.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const wallBoundedStreamLineParticle& p
)
{
    os  << static_cast<const wallBoundedParticle&>(p)
        << token::SPACE << p.lifeTime_
        << token::SPACE << p.sampledPositions_
        << token::SPACE << p.sampledScalars_
        << token::SPACE << p.sampledVectors_;

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
