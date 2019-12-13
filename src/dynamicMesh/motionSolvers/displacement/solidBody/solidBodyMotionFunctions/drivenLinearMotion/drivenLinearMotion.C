/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "drivenLinearMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "dimensionedVector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(drivenLinearMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        drivenLinearMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::drivenLinearMotion::drivenLinearMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
    CofGvelocity_(SBMFCoeffs.get<word>("CofGvelocity")),
    normal_(SBMFCoeffs.lookupOrDefault<vector>("normal", Zero)),
    CofGvel_
    (
        IOobject
        (
            CofGvelocity_,
            time_.timeName(),
            time_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        dimensionedVector(dimless, Zero)
    ),
    displacement_(Zero)
{
    read(SBMFCoeffs);
    if (mag(normal_) > SMALL)
    {
        normal_ /= (mag(normal_) + SMALL);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::drivenLinearMotion::~drivenLinearMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::drivenLinearMotion::transformation() const
{
    scalar deltaT = time_.deltaT().value();

    vector velocity = CofGvel_.value();

    // Take out normal component
    if (mag(normal_) > SMALL)
    {
        velocity = CofGvel_.value() - ((CofGvel_.value() & normal_)*normal_);
    }

    DebugInFunction << "Vel on plane  :" << velocity << endl;

    // Translation of centre of gravity with constant velocity
    //const vector displacement = velocity*t;
    displacement_ += velocity*deltaT;

    quaternion R(1);
    septernion TR(septernion(-displacement_)*R);

    DebugInFunction << "Time = " << time_.value()
                    << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::drivenLinearMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    return true;
}


// ************************************************************************* //
