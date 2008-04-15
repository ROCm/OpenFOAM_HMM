/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "engineTime.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from objectRegistry arguments
engineTime::engineTime
(
    const word& name,
    const fileName& rootPath,
    const fileName& caseName,
    const fileName& systemName,
    const fileName& constantName
)
:
    Time
    (
        name,
        rootPath,
        caseName,
        systemName,
        constantName
    ),
    engineGeometry_
    (
        IOobject
        (
            "engineGeometry",
            constant(),
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    conRodLength_(engineGeometry_.lookup("conRodLength")),
    bore_(engineGeometry_.lookup("bore")),
    stroke_(engineGeometry_.lookup("stroke")),
    clearance_(engineGeometry_.lookup("clearance")),
    rpm_(engineGeometry_.lookup("rpm"))
{
    value() = degToTime(value());

    startTime_ = degToTime(startTime_);
    endTime_ = degToTime(endTime_);

    deltaT_ = degToTime(deltaT_);
    deltaT0_ = deltaT_;

    if 
    (
        writeControl_ == wcRunTime
     || writeControl_ == wcAdjustableRunTime
    )
    {
        writeInterval_ = degToTime(writeInterval_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Read the controlDict and set all the parameters
bool engineTime::read()
{
    if (!Time::read())
    {
        return false;
    }
    else
    {
        deltaT_ = degToTime(deltaT_);
        endTime_ = degToTime(endTime_);

        if 
        (
            writeControl_ == wcRunTime
         || writeControl_ == wcAdjustableRunTime
        )
        {
            writeInterval_ = degToTime(writeInterval_);
        }

        return true;
    }
}


scalar engineTime::degToRad(const scalar deg) const
{
    return mathematicalConstant::pi*deg/180.0;
}


scalar engineTime::degToTime(const scalar theta) const
{
    return theta/(360.0*rpm_.value()/60.0);
}


scalar engineTime::timeToDeg(const scalar t) const
{
    return t*(360.0*rpm_.value()/60.0);
}


scalar engineTime::theta() const
{
    return timeToDeg(value());
}


// Return current crank-angle translated to a single revolution
// (value between -180 and 180 with 0 = top dead centre)
scalar engineTime::thetaRevolution() const
{
    scalar t = theta();

    while (t > 180.0)
    {
        t -= 360.0;
    }

    while (t < -180.0)
    {
        t += 360.0;
    }

    return t;
}


scalar engineTime::deltaTheta() const
{
    return timeToDeg(deltaT().value());
}


scalar engineTime::pistonPosition(const scalar theta) const
{
    return
    (
        conRodLength_.value()
      + stroke_.value()/2.0
      + clearance_.value()
    )
  - (
        stroke_.value()*::cos(degToRad(theta))/2.0
      + ::sqrt
        (
            sqr(conRodLength_.value())
            - sqr(stroke_.value()*::sin(degToRad(theta))/2.0)
        )
    );
}


dimensionedScalar engineTime::pistonPosition() const
{
    return dimensionedScalar
    (
        "pistonPosition",
        dimLength,
        pistonPosition(theta())
    );
}


dimensionedScalar engineTime::pistonDisplacement() const
{
    return dimensionedScalar
    (
        "pistonDisplacement",
        dimLength,
        pistonPosition(theta() - deltaTheta()) - pistonPosition().value()
    );
}


dimensionedScalar engineTime::pistonSpeed() const
{
    return dimensionedScalar
    (
        "pistonSpeed",
        dimLength/dimTime,
        pistonDisplacement().value()/(deltaT().value() + VSMALL)
    );
}


scalar engineTime::userTimeToTime(const scalar theta) const
{
    return degToTime(theta);
}


scalar engineTime::timeToUserTime(const scalar t) const
{
    return timeToDeg(t);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
