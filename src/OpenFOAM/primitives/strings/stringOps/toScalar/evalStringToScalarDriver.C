/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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

Description
    Ragel lexer interface for lemon grammar of a simple string to
    scalar evaluation

\*---------------------------------------------------------------------------*/

#include "evalStringToScalar.H"
#include "evalStringToScalarDriver.H"
#include "evalStringToScalarScanner.H"
#include "error.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parsing::evalStringToScalar::parseDriver::parseDriver()
:
    genericRagelLemonDriver(),
    value_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::parsing::evalStringToScalar::parseDriver::execute
(
    const std::string& s
)
{
    // scanner::debug = 1;

    scanner().process(s, *this);

    return value_;
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::stringOps::toScalar(const std::string& s)
{
    Foam::parsing::evalStringToScalar::parseDriver driver;

    scalar val = driver.execute(s);
    // val = driver.value();

    return val;
}


// ************************************************************************* //
