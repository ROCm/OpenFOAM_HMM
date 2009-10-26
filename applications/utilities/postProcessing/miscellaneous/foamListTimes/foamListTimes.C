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
    foamListTimes

Description
    List times using timeSelector

Usage

    - foamListTimes [OPTION]

    @param -processor \n
    List times from processor0 directory

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();  // -constant enabled
    argList::noBanner();
    argList::noParallel();
    argList::validOptions.insert("processor", "");

#   include "setRootCase.H"

    label nProcs = 0;

    // Create the processor databases
    PtrList<Time> databases(1);

    if (args.optionFound("processor"))
    {
        // determine the processor count directly
        while (isDir(args.path()/(word("processor") + name(nProcs))))
        {
            ++nProcs;
        }

        if (!nProcs)
        {
            FatalErrorIn(args.executable())
                << "No processor* directories found"
                << exit(FatalError);
        }

        // Create the processor databases
        databases.setSize(nProcs);

        forAll(databases, procI)
        {
            databases.set
            (
                procI,
                new Time
                (
                    Time::controlDictName,
                    args.rootPath(),
                    args.caseName()/fileName(word("processor") + name(procI))
                )
            );
        }
    }
    else
    {
        databases.set
        (
            0,
            new Time
            (
                Time::controlDictName,
                args.rootPath(),
                args.caseName()
            )
        );
    }


    // use the times list from the master processor
    // and select a subset based on the command-line options
    instantList timeDirs = timeSelector::select
    (
        databases[0].times(),
        args
    );

    forAll(timeDirs, timeI)
    {
        Info<< timeDirs[timeI].name() << endl;
    }

    return 0;
}


// ************************************************************************* //
