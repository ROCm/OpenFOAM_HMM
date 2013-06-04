/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Test-GAMGAgglomeration

Description
    Test application for GAMG agglomeration. Hardcoded to expect GAMG on p.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "GAMGAgglomeration.H"
#include "OFstream.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void checkConnectedAgglomeration
(
    const lduMesh& mesh,
    const labelUList& restrict,
    const label nCoarse
)
{
    if (mesh.lduAddr().size() != restrict.size())
    {
        FatalErrorIn
        (
            "checkConnectedAgglomeration(const lduMesh&, const labelList&)"
        )   << "nCells:" << mesh.lduAddr().size()
            << " agglom:" << restrict.size()
            << abort(FatalError);
    }

    // Seed (master) for every region
    labelList regionToMaster(nCoarse, -1);
    labelList master(mesh.lduAddr().size(), -1);
    forAll(restrict, cellI)
    {
        label region = restrict[cellI];
        if (regionToMaster[region] == -1)
        {
            // Set cell to be master for region
            //Pout<< "For region " << region
            //    << " allocating local master " << cellI
            //    << endl;
            regionToMaster[region] = cellI;
            master[cellI] = cellI;
        }
    }

    // Now loop and transport master through region
    const labelUList& lower = mesh.lduAddr().lowerAddr();
    const labelUList& upper = mesh.lduAddr().upperAddr();

    while (true)
    {
        label nChanged = 0;

        forAll(lower, faceI)
        {
            label own = lower[faceI];
            label nei = upper[faceI];

            if (restrict[own] == restrict[nei])
            {
                // Region-internal face

                if (master[own] != -1)
                {
                    if (master[nei] == -1)
                    {
                        master[nei] = master[own];
                        nChanged++;
                    }
                    else if (master[nei] != master[own])
                    {
                        FatalErrorIn("checkConnectedAgglomeration(..)")
                            << "problem" << abort(FatalError);
                    }
                }
                else if (master[nei] != -1)
                {
                    master[own] = master[nei];
                    nChanged++;
                }
            }
        }

        reduce(nChanged, sumOp<label>());

        if (nChanged == 0)
        {
            break;
        }
    }

    // Check that master is set for all cells
    boolList singleRegion(nCoarse, true);
    label nSet = nCoarse;
    forAll(master, cellI)
    {
        if (master[cellI] == -1)
        {
            label region = restrict[cellI];
            if (singleRegion[region] == true)
            {
                singleRegion[region] = false;
                nSet--;
            }
        }
    }

    label totalNCoarse = returnReduce(nCoarse, sumOp<label>());
    label totalNVisited = returnReduce(nSet, sumOp<label>());

    if (totalNVisited < totalNCoarse)
    {
        WarningIn("checkConnectedAgglomeration(..)")
            << "out of " << totalNCoarse
            << " agglomerated cells have " << totalNCoarse-totalNVisited
            << " cells that are not a single connected region" << endl;
    }
}


// Main program:

int main(int argc, char *argv[])
{
    argList::addBoolOption
    (
        "writeObj",
        "write obj files of agglomeration"
    );
    argList::addBoolOption
    (
        "normalise",
        "normalise agglomeration (0..1)"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    bool writeObj = args.optionFound("writeObj");
    bool normalise = args.optionFound("normalise");

    #include "createMesh.H"

    const fvSolution& sol = static_cast<const fvSolution&>(mesh);
    const dictionary& pDict = sol.subDict("solvers").subDict("p");

    const GAMGAgglomeration& agglom = GAMGAgglomeration::New
    (
        mesh,
        pDict
    );

    labelList cellToCoarse(identity(mesh.nCells()));
    labelListList coarseToCell(invertOneToMany(mesh.nCells(), cellToCoarse));

    runTime++;

    // Write initial agglomeration
    {
        volScalarField scalarAgglomeration
        (
            IOobject
            (
                "agglomeration",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("aggomeration", dimless, 0.0)
        );
        scalarField& fld = scalarAgglomeration.internalField();
        forAll(fld, cellI)
        {
            fld[cellI] = cellToCoarse[cellI];
        }
        fld /= max(fld);
        scalarAgglomeration.correctBoundaryConditions();
        scalarAgglomeration.write();

        Info<< "Writing initial cell distribution to "
            << runTime.timeName() << endl;
    }


    for (label level = 0; level < agglom.size(); level++)
    {
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        const labelList& addr = agglom.restrictAddressing(level);
        label coarseSize = max(addr)+1;

        Info<< "Level : " << level << endl
            << returnReduce(addr.size(), sumOp<label>()) << endl
            << "    current size      : "
            << returnReduce(addr.size(), sumOp<label>()) << endl
            << "    agglomerated size : "
            << returnReduce(coarseSize, sumOp<label>()) << endl;

        checkConnectedAgglomeration
        (
            agglom.meshLevel(level),
            addr,
            coarseSize
        );


        forAll(addr, fineI)
        {
            const labelList& cellLabels = coarseToCell[fineI];
            forAll(cellLabels, i)
            {
                cellToCoarse[cellLabels[i]] = addr[fineI];
            }
        }
        coarseToCell = invertOneToMany(coarseSize, cellToCoarse);

        // Write agglomeration
        {
            volScalarField scalarAgglomeration
            (
                IOobject
                (
                    "agglomeration",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("aggomeration", dimless, 0.0)
            );
            scalarField& fld = scalarAgglomeration.internalField();
            forAll(fld, cellI)
            {
                fld[cellI] = cellToCoarse[cellI];
            }
            if (normalise)
            {
                fld /= max(fld);
            }
            scalarAgglomeration.correctBoundaryConditions();
            scalarAgglomeration.write();
        }

        if (writeObj)
        {
            OFstream str(runTime.path()/runTime.timeName()/"aggomeration.obj");
            label vertI = 0;

            // Write all mesh cc
            forAll(mesh.cellCentres(), cellI)
            {
                meshTools::writeOBJ(str, mesh.cellCentres()[cellI]);
                vertI++;
            }

            // Determine coarse cc
            forAll(coarseToCell, coarseI)
            {
                const labelList& cellLabels = coarseToCell[coarseI];

                point coarseCc = average
                (
                    pointField(mesh.cellCentres(), cellLabels)
                );
                meshTools::writeOBJ(str, coarseCc);
                vertI++;

                forAll(cellLabels, i)
                {
                    label cellI = cellLabels[i];

                    str << "l " << cellI+1 << ' ' << vertI << nl;
                }
            }
        }

        Info<< endl;
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
