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
    checkDecomposePar

Group
    grpParallelUtilities

Description
    Check decomposition from kaffpa (KaHIP) output.
    foamToMetisGraph was likely used for producing the kaffpa input.

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "fvCFD.H"
#include "cpuTime.H"
#include "IFstream.H"
#include "regionProperties.H"
#include "decompositionInformation.H"
#include "decompositionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Check decomposition from kaffpa (KaHIP) output"
    );

    argList::noParallel();
    argList::noBanner();

    #include "addRegionOption.H"
    argList::addBoolOption
    (
        "allRegions",
        "operate on all regions in regionProperties"
    );
    argList::addBoolOption
    (
        "verbose",
        "more information about decomposition"
    );

    argList::addArgument("kaffpa-output-file");

    // Include explicit constant options, have zero from time range
    timeSelector::addOptions(true, false);

    #include "setRootCase.H"

    const fileName decompFile = args[1];

    const bool region     = args.optionFound("region");
    const bool allRegions = args.optionFound("allRegions");
    const bool verbose    = args.optionFound("verbose");

    // Set time from database
    #include "createTime.H"
    // Allow override of time
    instantList times = timeSelector::selectIfPresent(runTime, args);

    // Allow override of decomposeParDict location
    fileName decompDictFile;
    args.optionReadIfPresent("decomposeParDict", decompDictFile);

    wordList regionNames;
    wordList regionDirs;
    if (allRegions)
    {
        Info<< "Decomposing all regions in regionProperties" << nl << endl;
        regionProperties rp(runTime);
        forAllConstIters(rp, iter)
        {
            const wordList& regions = iter();
            forAll(regions, i)
            {
                if (!regionNames.found(regions[i]))
                {
                    regionNames.append(regions[i]);
                }
            }
        }
        regionDirs = regionNames;
    }
    else
    {
        word regionName;
        if (args.optionReadIfPresent("region", regionName))
        {
            regionNames = wordList(1, regionName);
            regionDirs = regionNames;
        }
        else
        {
            regionNames = wordList(1, fvMesh::defaultRegion);
            regionDirs = wordList(1, word::null);
        }
    }

    labelList cellToProc;

    forAll(regionNames, regioni)
    {
        const word& regionName = regionNames[regioni];
        const word& regionDir = regionDirs[regioni];

        Info<< "\n\nDecomposing mesh " << regionName << nl << endl;
        Info<< "Create mesh..." << flush;

        fvMesh mesh
        (
            IOobject
            (
                regionName,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        Info<< " nCells = " << mesh.nCells() << endl;

        // Expected format is a simple ASCII list
        cellToProc.setSize(mesh.nCells());
        {
            IFstream is(decompFile);

            forAll(cellToProc, celli)
            {
                cellToProc[celli] = readLabel(is);
            }
        }

        const label nDomains = max(cellToProc) + 1;

        CompactListList<label> cellCells;
        decompositionMethod::calcCellCells
        (
            mesh,
            identity(mesh.nCells()),
            mesh.nCells(),
            false,
            cellCells
        );

        decompositionInformation info
        (
            cellCells,
            cellToProc,
            nDomains
        );

        if (verbose)
        {
            info.printDetails(Info);
            Info<< nl;
        }
        info.printSummary(Info);
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
