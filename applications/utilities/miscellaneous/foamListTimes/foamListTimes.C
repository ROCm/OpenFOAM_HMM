/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
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
    List times using the timeSelector, or use to remove selected time
    directories.

Usage
    \b foamListTimes [OPTION]

    Options:
      - \par -processor
        List times from processor0/ directory

      - \par -rm
        Remove selected time directories

      - \par -verbose
        Report progress during removal

Note
    The OpenFOAM banner information is suppressed so that the output can be
    piped into another command.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "autoPtr.H"
#include "profiling.H"
#include "timeSelector.H"
#include "TimePaths.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "List times using the timeSelector,"
        " or use to remove selected time directories"
    );
    timeSelector::addOptions(true, true);  // constant(true), zero(true)
    argList::noBanner();
    argList::noParallel();
    argList::noJobInfo();
    argList::noFunctionObjects();  // Never use function objects
    argList::addBoolOption
    (
        "processor",
        "List times from processor0/ directory"
    );
    argList::addBoolOption
    (
        "rm",
        "Remove selected time directories"
    );
    argList::addBoolOption
    (
        "verbose",
        "Report progress of -rm option"
    );
    profiling::disable(); // Disable profiling (and its output)

    #include "setRootCase.H"

    const bool removeFiles(args.found("rm"));
    const bool verbose(args.found("verbose"));


    // Get times list from the master processor and subset based on
    // command-line options

    label nProcs = 0;
    autoPtr<TimePaths> timePaths;

    if (args.found("processor"))
    {
        // Determine the processor count
        nProcs = fileHandler().nProcs(args.path());

        if (!nProcs)
        {
            FatalErrorInFunction
                << "No processor* directories found"
                << exit(FatalError);
        }

        timePaths = autoPtr<TimePaths>::New
        (
            args.rootPath(),
            args.caseName()/"processor0"
        );
    }
    else
    {
        timePaths = autoPtr<TimePaths>::New
        (
            args.rootPath(),
            args.caseName()
        );
    }


    const instantList timeDirs(timeSelector::select(timePaths->times(), args));

    const label nTimes = timeDirs.size();

    if (removeFiles)
    {
        if (nProcs)
        {
            if (verbose)
            {
                Info<< "Removing " << nTimes
                    << " processor time directories" << endl;
            }

            forAllReverse(timeDirs, timei)
            {
                const word& timeName = timeDirs[timei].name();

                if (verbose)
                {
                    Info<< "    rm " << timeName
                        << " [" << (nTimes - timei) << '/' << nTimes << ']'
                        << endl;
                }

                fileName path(args.path()/"processors"/timeName);

                rmDir(path, true);

                for (label proci=0; proci<nProcs; ++proci)
                {
                    path =
                    (
                        args.path()
                      / ("processor" + Foam::name(proci))
                      / timeName
                    );

                    rmDir(path, true);
                }
            }
        }
        else
        {
            if (verbose)
            {
                Info<< "Removing " << nTimes
                    << " time directories" << endl;
            }

            forAllReverse(timeDirs, timei)
            {
                const word& timeName = timeDirs[timei].name();

                if (verbose)
                {
                    Info<< "    rm " << timeName
                        << " [" << (nTimes - timei) << '/' << nTimes << ']'
                        << endl;
                }

                rmDir(args.path()/timeName, true);
            }
        }
    }
    else
    {
        for (const instant& t : timeDirs)
        {
            Info<< t.name() << nl;
        }
        Info<< flush;
    }

    return 0;
}


// ************************************************************************* //
