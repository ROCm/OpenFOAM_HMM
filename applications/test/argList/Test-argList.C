/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOstreams.H"
#include "StringStream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    // argList::noFunctionObjects();
    argList::removeOption("case");

    argList::addOption("label", "value", "Test parsing of label");
    argList::addOption("scalar", "value", "Test parsing of scalar");

    // These are actually lies (never had -parseLabel, -parseScalar etc),
    // but good for testing...

    // Emits warning about it being old
    argList::addOptionCompat("label", {"parseLabel", 1612});

    // Specifying version=0 to use alias without any warnings
    argList::addOptionCompat("scalar", {"parseScalar", 0});

    // Fake a future option...
    argList::addOptionCompat("label", {"parse-label", 2112});

    argList args(argc, argv);

    Info<<"have: "
        <<args.optionCount({"label", "scalar"}) << " options" << nl;

    label ival;
    scalar sval;

    Info<< nl;

    Info<< "-label = " << flush;
    if (args.optionReadIfPresent("label", ival))
    {
        Info<< ival << endl;
    }
    else
    {
        Info<< "not specified" << endl;
    }

    Info<< "-scalar = " << flush;
    if (args.optionReadIfPresent("scalar", sval))
    {
        Info<< sval << endl;
    }
    else
    {
        Info<< "not specified" << endl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
