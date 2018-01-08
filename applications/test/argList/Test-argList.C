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
#include "Switch.H"
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
    argList::addOption("label",  "value", "Test parsing of label");
    argList::addOption("scalar", "value", "Test parsing of scalar");
    argList::addOption("string", "value", "Test string lookup");

    // These are actually lies (never had -parseLabel, -parseScalar etc),
    // but good for testing...

    // Emits warning about it being old
    argList::addOptionCompat("label", {"parseLabel", 1612});

    // Specifying version=0 to use alias without any warnings
    argList::addOptionCompat("scalar", {"parseScalar", 0});

    // Fake a future option...
    argList::addOptionCompat("label", {"parse-label", 2112});

    argList::addArgument("label");
    argList::addArgument("...");
    argList::addArgument("label");
    argList::nonMandatoryArgs();

    argList args(argc, argv, false, true);

    Info<<"have: "
        <<args.count({"label", "scalar"}) << " options" << nl;

    label ival;
    scalar sval;

    Info<< nl;

    Info<< "-label = " << flush;
    if (args.readIfPresent("label", ival))
    {
        Info<< ival << nl;
    }
    else
    {
        Info<< "not specified" << nl;
    }

    Info<< "-scalar = " << flush;
    if (args.readIfPresent("scalar", sval))
    {
        Info<< sval << nl;
    }
    else
    {
        Info<< "not specified" << nl;
    }


    // Using direct reading
    Info<< nl;
    if (args.found("label"))
    {
        Info<< "-label = " << args.opt<label>("label")
            << " or " << args.optionRead<label>("label")  // old-compat
            << " or " << readLabel(args["label"])         // with function
            << nl;
    }

    if (args.found("scalar"))
    {
        Info<< "-scalar = " << args.opt<scalar>("scalar")
            << " or " << args.optionRead<scalar>("scalar") // old-compat
            << " or " << readScalar(args["scalar"])        // with function
            << nl;
    }

    if (args.found("string"))
    {
        Info<< "-string = " << args.opt("string")
            << " or " << args.optionRead<scalar>("string")  // old-compat
            << nl;
    }


    // Arg reading
    Info<< nl;
    for (label argi=1; argi < args.size(); ++argi)
    {
        Info<< "arg[" << argi << "] = " << args.read<string>(argi)
            << " or " << args.argRead<label>(argi)  // old-compat
            << nl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
