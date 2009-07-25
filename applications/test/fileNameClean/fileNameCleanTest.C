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

Application
    fileNameCleanTest

Description


\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fileName.H"
#include "SubList.H"
#include "IOobject.H"
#include "IOstreams.H"
#include "OSspecific.H"


using namespace Foam;

void printCleaning(fileName& pathName)
{
    Info<< "fileName = " << pathName << nl
        << "  path() = " << pathName.path() << nl
        << "  name() = " << pathName.name() << nl
        << "  joined = " << pathName.path()/pathName.name() << nl << nl;

    pathName.clean();

    Info<< "cleaned  = " << pathName << nl
        << "  path() = " << pathName.path() << nl
        << "  name() = " << pathName.name() << nl
        << "  joined = " << pathName.path()/pathName.name() << nl << nl;

    IOobject::writeDivider(Info);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::validArgs.insert("fileName .. fileNameN");

    argList args(argc, argv, false, true);

    if (args.additionalArgs().empty())
    {
        args.printUsage();
    }

    if (args.optionFound("case"))
    {
        fileName pathName = args.option("case");
        Info<< nl
            << "-case" << nl
            << "path = " << args.path() << nl
            << "root = " << args.rootPath() << nl
            << "case = " << args.caseName() << nl
            << "FOAM_CASE=" << getEnv("FOAM_CASE") << nl
            << "FOAM_CASENAME=" << getEnv("FOAM_CASENAME") << nl
            << endl;

        printCleaning(pathName);
    }

    forAll(args.additionalArgs(), argI)
    {
        fileName pathName = args.additionalArgs()[argI];
        printCleaning(pathName);
    }

    Info<< "\nEnd" << endl;

    return 0;
}


// ************************************************************************* //
