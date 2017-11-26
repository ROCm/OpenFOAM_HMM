/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

Application
    foamListTimes

Group
    grpPostProcessingUtilities

Description
    List times using timeSelector.
    To simplify parsing of the output, the normal banner information
    is suppressed.

Usage
    \b foamListTimes [OPTION]

    Options:
      - \par -rm
        Remove selected time directories

      - \par -processor
        List times from processor0/ directory

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "profiling.H"
#include "timeSelector.H"
#include "Time.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote("List times using timeSelector");

    timeSelector::addOptions(true, true);
    argList::noBanner();
    argList::noParallel();
    argList::noJobInfo();
    argList::noFunctionObjects();
    argList::addBoolOption
    (
        "processor",
        "list times from processor0/ directory"
    );
    argList::addBoolOption
    (
        "rm",
        "remove selected time directories"
    );
    profiling::disable(); // Disable profiling (and its output)

    #include "setRootCase.H"

    // Get times list from the master processor and subset based on
    // command-line options

    label nProcs = 0;
    instantList timeDirs;

    if (args.optionFound("processor"))
    {
        // Determine the processor count
        nProcs = fileHandler().nProcs(args.path());

        if (!nProcs)
        {
            FatalErrorInFunction
                << "No processor* directories found"
                << exit(FatalError);
        }

        timeDirs = timeSelector::select
        (
            Time
            (
                Time::controlDictName,
                args.rootPath(),
                args.caseName()/"processor0"
            ).times(),
            args
        );
    }
    else
    {
        timeDirs = timeSelector::select
        (
            Time
            (
                Time::controlDictName,
                args.rootPath(),
                args.caseName()
            ).times(),
            args
        );
    }


    if (args.optionFound("rm"))
    {
        if (nProcs)
        {
            // Info<< "Remove " << timeDirs.size()
            //     << " processor time directories" << nl;

            forAllReverse(timeDirs, timei)
            {
                fileName path
                (
                    args.path()
                  / "processors"
                  / timeDirs[timei].name()
                );

                rmDir(path, true);

                for (label proci=0; proci<nProcs; ++proci)
                {
                    path =
                    (
                        args.path()
                      / (word("processor") + name(proci))
                      / timeDirs[timei].name()
                    );

                    rmDir(path, true);
                }
            }
        }
        else
        {
            // Info<< "Remove " << timeDirs.size()
            //     << " time directories" << nl;

            forAllReverse(timeDirs, timei)
            {
                rmDir(args.path()/timeDirs[timei].name(), true);
            }
        }
    }
    else
    {
        forAll(timeDirs, timei)
        {
            Info<< timeDirs[timei].name() << nl;
        }
        Info<< flush;
    }


    return 0;
}


// ************************************************************************* //
