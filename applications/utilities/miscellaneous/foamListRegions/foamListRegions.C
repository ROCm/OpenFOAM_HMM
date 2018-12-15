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

Application
    foamListRegions

Group
    grpPostProcessingUtilities

Description
    List regions from constant/regionProperties.

Usage
    \b foamListRegions [OPTION]

Note
    The OpenFOAM banner information is suppressed so that the output can be
    piped into another command.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "regionProperties.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "List regions from constant/regionProperties"
    );

    argList::noBanner();
    argList::noParallel();
    argList::noJobInfo();
    argList::noFunctionObjects();  // Never use function objects
    // No profiling since there is no time loop

    // Arguments are optional (non-mandatory)
    argList::noMandatoryArgs();
    argList::addArgument("regionType ... regionType");

    #include "setRootCase.H"

    // As per "createTime.H", but quieter.
    Time runTime(Time::controlDictName, args);

    regionProperties rp(runTime);

    // We now handle checking args and general sanity etc.
    wordList regionTypes;

    if (args.size() > 1)
    {
        regionTypes.setSize(args.size()-1);

        label nTypes = 0;
        for (label argi = 1; argi < args.size(); ++argi)
        {
            regionTypes[nTypes] = args[argi];

            if (rp.found(regionTypes[nTypes]))
            {
                ++nTypes;
            }
            else
            {
                std::cerr<< "No region-type: " << regionTypes[nTypes] << nl;
            }
        }

        regionTypes.setSize(nTypes);
    }
    else
    {
        regionTypes = rp.sortedToc();
    }


    for (const word& regionType : regionTypes)
    {
        const wordList& regionNames = rp[regionType];

        for (const word& regionName : regionNames)
        {
            Info<< regionName << nl;
        }
    }

    return 0;
}


// ************************************************************************* //
