/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

    #include "addAllRegionOptions.H"

    argList::addVerboseOption
    (
        "more information about decomposition"
    );

    argList::addArgument("kaffpa-output-file");

    // Include explicit constant options, have zero from time range
    timeSelector::addOptions(true, false);

    #include "setRootCase.H"

    const auto decompFile = args.get<fileName>(1);

    // Set time from database
    #include "createTime.H"

    // Allow override of time
    instantList times = timeSelector::selectIfPresent(runTime, args);

    // Allow override of decomposeParDict location
    const fileName decompDictFile =
        args.getOrDefault<fileName>("decomposeParDict", "");

    // Get region names
    #include "getAllRegionOptions.H"

    labelList cellToProc;

    forAll(regionNames, regioni)
    {
        const word& regionName = regionNames[regioni];
        // const word& regionDir =
        // (
        //     regionName != polyMesh::defaultRegion
        //   ? regionName
        //   : word::null
        // );

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

        if (args.verbose())
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
