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

\*----------------------------------------------------------------------------*/

#include "streamLineParticle.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(streamLineParticle, 0);
    defineParticleTypeNameAndDebug(streamLineParticle, 0);
    defineTemplateTypeNameAndDebugWithName
    (
        IOField<vectorField>,
        "vectorFieldField",
        0
    );
};

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::vector Foam::streamLineParticle::interpolateFields
(
    const streamLineParticle::trackData& td,
    const point& position,
    const label cellI
)
{
    if (cellI == -1)
    {
        FatalErrorIn("streamLineParticle::interpolateFields(..)")
            << "Cell:" << cellI << abort(FatalError);
    }

    sampledScalars_.setSize(td.vsInterp_.size());
    forAll(td.vsInterp_, scalarI)
    {
        sampledScalars_[scalarI].append
        (
            td.vsInterp_[scalarI].interpolate
            (
                position,
                cellI
            )
        );
    }

    sampledVectors_.setSize(td.vvInterp_.size());
    forAll(td.vvInterp_, vectorI)
    {
        sampledVectors_[vectorI].append
        (
            td.vvInterp_[vectorI].interpolate
            (
                position,
                cellI
            )
        );
    }

    const DynamicList<vector>& U = sampledVectors_[td.UIndex_];

    return U[U.size()-1];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
Foam::streamLineParticle::streamLineParticle
(
    const Cloud<streamLineParticle>& c,
    const vector& position,
    const label celli,
    const label lifeTime
)
:
    Particle<streamLineParticle>(c, position, celli),
    lifeTime_(lifeTime)
{}


//- Construct from Istream
Foam::streamLineParticle::streamLineParticle
(
    const Cloud<streamLineParticle>& c,
    Istream& is,
    bool readFields
)
:
    Particle<streamLineParticle>(c, is, readFields)
{
    if (readFields)
    {
        //if (is.format() == IOstream::ASCII)
        List<scalarList> sampledScalars;
        List<vectorList> sampledVectors;

        is  >> lifeTime_ >> sampledPositions_ >> sampledScalars
            >> sampledVectors;

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

    // Check state of Istream
    is.check
    (
        "streamLineParticle::streamLineParticle"
        "(const Cloud<streamLineParticle>&, Istream&, bool)"
    );
}


// Construct copy
Foam::streamLineParticle::streamLineParticle
(
    const streamLineParticle& c
)
:
    Particle<streamLineParticle>(c),
    lifeTime_(c.lifeTime_),
    sampledPositions_(c.sampledPositions_),
    sampledScalars_(c.sampledScalars_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::streamLineParticle::move(streamLineParticle::trackData& td)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    scalar deltaT = GREAT;  //cloud().pMesh().time().deltaT().value();
    scalar tEnd = (1.0 - stepFraction())*deltaT;
    scalar dtMax = tEnd;

    while
    (
        td.keepParticle
    && !td.switchProcessor
    && lifeTime_ > 0
    && tEnd > ROOTVSMALL
    )
    {
        // set the lagrangian time-step
        scalar dt = min(dtMax, tEnd);

        // Store current position and sampled velocity.
        --lifeTime_;
        sampledPositions_.append(position());
        vector U = interpolateFields(td, position(), cell());

        if (!td.trackForward_)
        {
            U = -U;
        }

        dt *= trackToFace(position()+dt*U, td);

        tEnd -= dt;
        stepFraction() = 1.0 - tEnd/deltaT;
    }

    if (!td.keepParticle || lifeTime_ == 0)
    {
        if (lifeTime_ == 0)
        {
            if (debug)
            {
                Pout<< "streamLineParticle : Removing stagnant particle:"
                    << static_cast<Particle<streamLineParticle> >(*this)
                    << " sampled positions:" << sampledPositions_.size()
                    << endl;
            }
            td.keepParticle = false;
        }
        else
        {
            // Normal exit. Store last position and fields
            sampledPositions_.append(position());
            interpolateFields(td, position(), cell());

            if (debug)
            {
                Pout<< "streamLineParticle : Removing particle:"
                    << static_cast<Particle<streamLineParticle> >(*this)
                    << " sampled positions:" << sampledPositions_.size()
                    << endl;
            }
        }

        // Transfer particle data into trackData.
        //td.allPositions_.append(sampledPositions_);
        td.allPositions_.append(vectorList());
        vectorList& top = td.allPositions_[td.allPositions_.size()-1];
        top.transfer(sampledPositions_);

        forAll(sampledScalars_, i)
        {
            //td.allScalars_[i].append(sampledScalars_[i]);
            td.allScalars_[i].append(scalarList());
            scalarList& top = td.allScalars_[i][td.allScalars_[i].size()-1];
            top.transfer(sampledScalars_[i]);
        }
        forAll(sampledVectors_, i)
        {
            //td.allVectors_[i].append(sampledVectors_[i]);
            td.allVectors_[i].append(vectorList());
            vectorList& top = td.allVectors_[i][td.allVectors_[i].size()-1];
            top.transfer(sampledVectors_[i]);
        }
    }

    return td.keepParticle;
}


bool Foam::streamLineParticle::hitPatch
(
    const polyPatch&,
    streamLineParticle::trackData& td,
    const label patchI
)
{
    // Disable generic patch interaction
    return false;
}


bool Foam::streamLineParticle::hitPatch
(
    const polyPatch&,
    int&,
    const label
)
{
    // Disable generic patch interaction
    return false;
}


void Foam::streamLineParticle::hitWedgePatch
(
    const wedgePolyPatch& pp,
    streamLineParticle::trackData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamLineParticle::hitWedgePatch
(
    const wedgePolyPatch&,
    int&
)
{}


void Foam::streamLineParticle::hitSymmetryPatch
(
    const symmetryPolyPatch& pp,
    streamLineParticle::trackData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamLineParticle::hitSymmetryPatch
(
    const symmetryPolyPatch&,
    int&
)
{}


void Foam::streamLineParticle::hitCyclicPatch
(
    const cyclicPolyPatch& pp,
    streamLineParticle::trackData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamLineParticle::hitCyclicPatch
(
    const cyclicPolyPatch&,
    int&
)
{}


void Foam::streamLineParticle::hitProcessorPatch
(
    const processorPolyPatch&,
    streamLineParticle::trackData& td
)
{
    // Switch particle
    td.switchProcessor = true;
}


void Foam::streamLineParticle::hitProcessorPatch
(
    const processorPolyPatch&,
    int&
)
{}


void Foam::streamLineParticle::hitWallPatch
(
    const wallPolyPatch& wpp,
    streamLineParticle::trackData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamLineParticle::hitWallPatch
(
    const wallPolyPatch& wpp,
    int&
)
{}


void Foam::streamLineParticle::hitPatch
(
    const polyPatch& wpp,
    streamLineParticle::trackData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamLineParticle::hitPatch
(
    const polyPatch& wpp,
    int&
)
{}


void Foam::streamLineParticle::readFields(Cloud<streamLineParticle>& c)
{
    if (!c.size())
    {
        return;
    }

    IOField<label> lifeTime
    (
        c.fieldIOobject("lifeTime", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, lifeTime);

    IOField<pointField> sampledPositions
    (
        c.fieldIOobject("sampledPositions", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, sampledPositions);

//    IOField<pointField> sampleVelocity
//    (
//        c.fieldIOobject("sampleVelocity", IOobject::MUST_READ)
//    );
//    c.checkFieldIOobject(c, sampleVelocity);

    label i = 0;
    forAllIter(Cloud<streamLineParticle>, c, iter)
    {
        iter().lifeTime_ = lifeTime[i];
        iter().sampledPositions_.transfer(sampledPositions[i]);
//        iter().sampleVelocity_.transfer(sampleVelocity[i]);
        i++;
    }
}


void Foam::streamLineParticle::writeFields(const Cloud<streamLineParticle>& c)
{
    Particle<streamLineParticle>::writeFields(c);

    label np =  c.size();

    IOField<label> lifeTime
    (
        c.fieldIOobject("lifeTime", IOobject::NO_READ),
        np
    );
    IOField<pointField> sampledPositions
    (
        c.fieldIOobject("sampledPositions", IOobject::NO_READ),
        np
    );
//    IOField<pointField> sampleVelocity
//    (
//        c.fieldIOobject("sampleVelocity", IOobject::NO_READ),
//        np
//    );

    label i = 0;
    forAllConstIter(Cloud<streamLineParticle>, c, iter)
    {
        lifeTime[i] = iter().lifeTime_;
        sampledPositions[i] = iter().sampledPositions_;
//        sampleVelocity[i] = iter().sampleVelocity_;
        i++;
    }

    lifeTime.write();
    sampledPositions.write();
//    sampleVelocity.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const streamLineParticle& p)
{
    os  << static_cast<const Particle<streamLineParticle>&>(p)
        << token::SPACE << p.lifeTime_
        << token::SPACE << p.sampledPositions_
        << token::SPACE << p.sampledScalars_
        << token::SPACE << p.sampledVectors_;

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const streamLineParticle&)");

    return os;
}


// ************************************************************************* //
