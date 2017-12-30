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
    decomposePar

Group
    grpParallelUtilities

Description
    Automatically decomposes a mesh and fields of a case for parallel
    execution of OpenFOAM.

Usage
    \b decomposePar [OPTION]

    Options:
      - \par -region \<regionName\>
        Decompose named region. Does not check for existence of processor*.

      - \par -allRegions
        Decompose all regions in regionProperties. Does not check for
        existence of processor*.

      - \par -constant

      - \par -time xxx:yyy
        Override controlDict settings and decompose selected times. Does not
        re-decompose the mesh i.e. does not handle moving mesh or changing
        mesh cases.

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
        "decompose a mesh and fields of a case for parallel execution"
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
        "operate on all regions in regionProperties"
    );
    argList::addBoolOption
    (
        "verbose",
        "more information about decomposition"
    );
    argList::addOption
    (
        "domains",
        "N"
        "override numberOfSubdomains"
    );

    argList::addOption
    (
        "method",
        "name"
        "override method"
    );

    // Include explicit constant options, have zero from time range
    timeSelector::addOptions(true, false);

    #include "setRootCase.H"

    const bool region     = args.optionFound("region");
    const bool allRegions = args.optionFound("allRegions");
    const bool verbose    = args.optionFound("verbose");

    const label numSubdomains =
        args.optionLookupOrDefault<label>("domains", 0);

    const word methodName =
        args.optionLookupOrDefault<word>("method", word::null);


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
