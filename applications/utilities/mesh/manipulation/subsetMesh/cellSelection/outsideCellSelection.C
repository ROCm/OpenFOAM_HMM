/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "outsideCellSelection.H"
#include "addToRunTimeSelectionTable.H"
#include "faceSet.H"
#include "polyMesh.H"
#include "motionSmoother.H"
#include "regionSplit.H"
#include "syncTools.H"
#include "zeroGradientFvPatchFields.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cellSelections
{
    defineTypeNameAndDebug(outsideCellSelection, 0);
    addToRunTimeSelectionTable(cellSelection, outsideCellSelection, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::cellSelections::outsideCellSelection::generateField
(
    const word& name,
    const boolList& lst
) const
{
    const fvMesh& mesh = dynamic_cast<const fvMesh&>(mesh_);

    tmp<volScalarField> tfld
    (
        new volScalarField
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(name, dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    scalarField& fld = tfld().internalField();

    forAll(fld, celli)
    {
       fld[celli] = 1.0*lst[celli];
    }
    tfld().correctBoundaryConditions();

    return tfld;
}


void Foam::cellSelections::outsideCellSelection::markRegionFaces
(
    const boolList& selectedCell,
    boolList& regionFace
) const
{
    // Internal faces
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();
    forAll(faceNeighbour, faceI)
    {
        if
        (
            selectedCell[faceOwner[faceI]]
         != selectedCell[faceNeighbour[faceI]]
        )
        {
            regionFace[faceI] = true;
        }
    }

    // Swap neighbour selectedCell state
    boolList nbrSelected;
    syncTools::swapBoundaryCellList(mesh_, selectedCell, nbrSelected);

    // Boundary faces
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];
        const labelUList& faceCells = pp.faceCells();
        forAll(faceCells, i)
        {
            label faceI = pp.start()+i;
            label bFaceI = faceI-mesh_.nInternalFaces();
            if
            (
                selectedCell[faceCells[i]]
             != selectedCell[nbrSelected[bFaceI]]
            )
            {
                regionFace[faceI] = true;
            }
        }
    }
}


Foam::boolList Foam::cellSelections::outsideCellSelection::findRegions
(
    const bool verbose,
    const regionSplit& cellRegion
) const
{
    boolList keepRegion(cellRegion.nRegions(), false);

    forAll(locationsInMesh_, i)
    {
        // Find the region containing the insidePoint

        label cellI = mesh_.findCell(locationsInMesh_[i]);

        label keepRegionI = -1;
        label keepProcI = -1;
        if (cellI != -1)
        {
            keepRegionI = cellRegion[cellI];
            keepProcI = Pstream::myProcNo();
        }
        reduce(keepRegionI, maxOp<label>());
        keepRegion[keepRegionI] = true;

        reduce(keepProcI, maxOp<label>());

        if (keepProcI == -1)
        {
            FatalErrorIn
            (
                "outsideCellSelection::findRegions"
                "(const bool, const regionSplit&)"
            )   << "Did not find " << locationsInMesh_[i]
                << " in mesh." << " Mesh bounds are " << mesh_.bounds()
                << exit(FatalError);
        }

        if (verbose)
        {
            Info<< "Found location " << locationsInMesh_[i]
                << " in cell " << cellI << " on processor " << keepProcI
                << " in global region " << keepRegionI
                << " out of " << cellRegion.nRegions() << " regions." << endl;
        }
    }

    return keepRegion;
}


void Foam::cellSelections::outsideCellSelection::unselectOutsideRegions
(
    boolList& selectedCell
) const
{
    // Determine faces on the edge of selectedCell
    boolList blockedFace(mesh_.nFaces(), false);
    markRegionFaces(selectedCell, blockedFace);

    // Determine regions
    regionSplit cellRegion(mesh_, blockedFace);

    // Determine regions containing locationsInMesh_
    boolList keepRegion(findRegions(true, cellRegion));

    // Go back to bool per cell
    forAll(cellRegion, cellI)
    {
        if (!keepRegion[cellRegion[cellI]])
        {
            selectedCell[cellI] = false;
        }
    }
}


void Foam::cellSelections::outsideCellSelection::shrinkRegions
(
    boolList& selectedCell
) const
{
    // Select points on unselected cells and boundary
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    boolList boundaryPoint(mesh_.nPoints(), false);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];

        if (!pp.coupled() && !isA<emptyPolyPatch>(pp))
        {
            forAll(pp, i)
            {
                const face& f = pp[i];
                forAll(f, fp)
                {
                    boundaryPoint[f[fp]] = true;
                }
            }
        }
    }

    forAll(selectedCell, cellI)
    {
        if (!selectedCell[cellI])
        {
            const labelList& cPoints = mesh_.cellPoints(cellI);
            forAll(cPoints, i)
            {
                boundaryPoint[cPoints[i]] = true;
            }
        }
    }

    syncTools::syncPointList(mesh_, boundaryPoint, orEqOp<bool>(), false);


    // Select all cells using these points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    label nChanged = 0;
    forAll(boundaryPoint, pointI)
    {
        if (boundaryPoint[pointI])
        {
            const labelList& pCells = mesh_.pointCells(pointI);
            forAll(pCells, i)
            {
                label cellI = pCells[i];
                if (selectedCell[cellI])
                {
                    selectedCell[cellI] = false;
                    nChanged++;
                }
            }
        }
    }
}


void Foam::cellSelections::outsideCellSelection::erode
(
    boolList& selectedCell
) const
{
    //Info<< "Entering shrinkRegions:" << count(selectedCell) << endl;
    //generateField("selectedCell_before", selectedCell)().write();

    // Now erode and see which regions get disconnected
    boolList shrunkSelectedCell(selectedCell);

    for (label iter = 0; iter < nErode_; iter++)
    {
        shrinkRegions(shrunkSelectedCell);
    }

    //Info<< "After shrinking:" << count(shrunkSelectedCell) << endl;
    //generateField("shrunkSelectedCell", shrunkSelectedCell)().write();



    // Determine faces on the edge of shrunkSelectedCell
    boolList blockedFace(mesh_.nFaces(), false);
    markRegionFaces(shrunkSelectedCell, blockedFace);

    // Find disconnected regions
    regionSplit cellRegion(mesh_, blockedFace);

    // Determine regions containing insidePoints
    boolList keepRegion(findRegions(true, cellRegion));


    // Extract cells in regions that are not to be kept.
    boolList removeCell(mesh_.nCells(), false);
    forAll(cellRegion, cellI)
    {
        if (shrunkSelectedCell[cellI] && !keepRegion[cellRegion[cellI]])
        {
            removeCell[cellI] = true;
        }
    }

    //Info<< "removeCell before:" << count(removeCell) << endl;
    //generateField("removeCell_before", removeCell)().write();



    // Grow removeCell
    for (label iter = 0; iter < nErode_; iter++)
    {
        // Grow selected cell in regions that are not for keeping
        boolList boundaryPoint(mesh_.nPoints(), false);
        forAll(removeCell, cellI)
        {
            if (removeCell[cellI])
            {
                const labelList& cPoints = mesh_.cellPoints(cellI);
                forAll(cPoints, i)
                {
                    boundaryPoint[cPoints[i]] = true;
                }
            }
        }
        syncTools::syncPointList(mesh_, boundaryPoint, orEqOp<bool>(), false);

        // Select all cells using these points

        label nChanged = 0;
        forAll(boundaryPoint, pointI)
        {
            if (boundaryPoint[pointI])
            {
                const labelList& pCells = mesh_.pointCells(pointI);
                forAll(pCells, i)
                {
                    label cellI = pCells[i];
                    if (!removeCell[cellI])
                    {
                        removeCell[cellI] = true;
                        nChanged++;
                    }
                }
            }
        }
    }

    //Info<< "removeCell after:" << count(removeCell) << endl;
    //generateField("removeCell_after", removeCell)().write();


    // Unmark removeCell
    forAll(removeCell, cellI)
    {
        if (removeCell[cellI])
        {
            selectedCell[cellI] = false;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellSelections::outsideCellSelection::outsideCellSelection
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    cellSelection(name, mesh, dict),
    locationsInMesh_(dict.lookup("locationsInMesh")),
    nErode_(readLabel(dict.lookup("nErodeLayers")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellSelections::outsideCellSelection::~outsideCellSelection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cellSelections::outsideCellSelection::select
(
    boolList& selectedCell
) const
{
    // Unselect all disconnected regions
    unselectOutsideRegions(selectedCell);

    if (nErode_ > 0)
    {
        erode(selectedCell);
    }
}


// ************************************************************************* //
