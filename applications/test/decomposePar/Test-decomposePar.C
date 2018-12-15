/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenCFD Ltd.
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
    Test-decomposePar

Group
    grpParallelUtilities

Description
    Like decomposePar -dry-run, but with additional options

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "fvCFD.H"
#include "cpuTime.H"
#include "IOobjectList.H"
#include "regionProperties.H"
#include "decompositionInformation.H"
#include "decompositionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Special-purpose version of decomposePar with additional"
        " -domain and -method options."
        " The '-dry-run' and '-cellDist' are implicit.\n"
        "NB: The -domain/-method overrides may not work very well with regions"
    );

    argList::noParallel();
    argList::addOption
    (
        "decomposeParDict",
        "file",
        "read decomposePar dictionary from specified location"
    );
    #include "addRegionOption.H"
    argList::addBoolOption
    (
        "allRegions",
        "Operate on all regions in regionProperties"
    );
    argList::addBoolOption
    (
        "verbose",
        "Additional verbosity"
    );
    argList::addOption
    (
        "domains",
        "N",
        "Override numberOfSubdomains"
    );

    argList::addOption
    (
        "method",
        "name",
        "Override decomposition method"
    );


    // These are implicit so just ignore them
    argList::ignoreOptionCompat({"dry-run", 0}, false);
    argList::ignoreOptionCompat({"cellDist", 0}, false);

    // Include explicit constant options, have zero from time range
    timeSelector::addOptions(true, false);

    #include "setRootCase.H"

    const bool optRegion  = args.found("region");
    const bool allRegions = args.found("allRegions");
    const bool verbose    = args.found("verbose");

    const label numSubdomains = args.opt<label>("domains", 0);
    const word methodName = args.opt<word>("method", word::null);

    // Set time from database
    #include "createTime.H"
    // Allow override of time
    instantList times = timeSelector::selectIfPresent(runTime, args);

    // Allow override of decomposeParDict location
    const fileName decompDictFile = args.opt<fileName>("decomposeParDict", "");

    // Get all region names
    wordList regionNames;
    if (allRegions)
    {
        regionNames = regionProperties(runTime).names();

        Info<< "Decomposing all regions in regionProperties" << nl
            << "    " << flatOutput(regionNames) << nl << endl;
    }
    else
    {
        regionNames.resize(1);
        regionNames.first() = args.opt<word>("region", fvMesh::defaultRegion);
    }

    forAll(regionNames, regioni)
    {
        const word& regionName = regionNames[regioni];
        const word& regionDir =
            (regionName == fvMesh::defaultRegion ? word::null : regionName);

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

        Info<< "\nCalculating distribution of cells" << endl;
        cpuTime decompositionTime;

        const decompositionModel& model = decompositionModel::New
        (
            mesh,
            decompDictFile
        );

        // Allow command-line override for quick testing

        dictionary& modelDict = const_cast<decompositionModel&>(model);

        if (numSubdomains)
        {
            modelDict.add
            (
                word("numberOfSubdomains"),
                numSubdomains,
                true
            );
        }

        if (!methodName.empty())
        {
            modelDict.add
            (
                word("method"),
                methodName,
                true
            );
        }

        scalarField cellWeights;
        word weightName;
        if (model.readIfPresent("weightField", weightName))
        {
            volScalarField weights
            (
                IOobject
                (
                    weightName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );
            cellWeights = weights.primitiveField();
        }

        decompositionMethod& method = model.decomposer();

        CompactListList<label> cellCells;
        decompositionMethod::calcCellCells
        (
            mesh,
            identity(mesh.nCells()),
            mesh.nCells(),
            false,
            cellCells
        );

        labelList cellToProc = method.decompose(mesh, cellWeights);

        Info<< "\nFinished decomposition into "
            << method.nDomains() << " domains in "
            << decompositionTime.elapsedCpuTime()
            << " s" << nl << endl;

        decompositionInformation info
        (
            cellCells,
            cellToProc,
            method.nDomains()
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
