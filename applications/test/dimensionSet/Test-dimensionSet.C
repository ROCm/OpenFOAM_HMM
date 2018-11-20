/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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
    Print values of predefined dimensionSets, and some other tests

\*---------------------------------------------------------------------------*/

#include "dimensionSet.H"
#include "IOmanip.H"

using namespace Foam;

// Macros to stringify macro contents.
#define STRINGIFY(content)      #content
#define STRING_QUOTE(input)     STRINGIFY(input)

#define PRINT_DIMSET(arg)       \
    Info<< STRING_QUOTE(arg) << "  " << arg << nl


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    Info<< "Predefined dimensionSets" << nl << nl;

    PRINT_DIMSET(dimless);
    PRINT_DIMSET(dimMass);
    PRINT_DIMSET(dimLength);
    PRINT_DIMSET(dimTime);
    PRINT_DIMSET(dimTemperature);
    PRINT_DIMSET(dimMoles);
    PRINT_DIMSET(dimCurrent);
    PRINT_DIMSET(dimLuminousIntensity);

    PRINT_DIMSET(dimArea);
    PRINT_DIMSET(dimVolume);
    PRINT_DIMSET(dimVol);

    PRINT_DIMSET(dimVelocity);
    PRINT_DIMSET(dimAcceleration);

    PRINT_DIMSET(dimDensity);
    PRINT_DIMSET(dimForce);
    PRINT_DIMSET(dimEnergy);
    PRINT_DIMSET(dimPower);

    PRINT_DIMSET(dimPressure);
    PRINT_DIMSET(dimCompressibility);
    PRINT_DIMSET(dimGasConstant);
    PRINT_DIMSET(dimSpecificHeatCapacity);
    PRINT_DIMSET(dimViscosity);
    PRINT_DIMSET(dimDynamicViscosity);


    Info<< nl << "Operators" << nl << nl;

    // Operators
    {
        Info<< "dimLength/dimTime  " << (dimLength/dimTime) << nl;
        Info<< "1/dimTime  " << (dimless/dimTime) << nl;
        Info<< "-dimTime  " << (-dimTime) << nl;
        Info<< "~dimTime  " << (~dimTime) << nl;

        Info<< "sqr(dimVelocity) " << sqr(dimVelocity) << nl;
        Info<< "pow2(dimVelocity) " << pow2(dimVelocity) << nl;
        Info<< "sqrt(dimVelocity) " << sqrt(dimVelocity) << nl;
        Info<< "pow025(dimVelocity) " << pow025(dimVelocity) << nl;
        Info<< "sqrt(sqrt(dimVelocity)) " << sqrt(sqrt(dimVelocity)) << nl;
    }


    Info<< nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
