/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "boundaryFirstRenumber.H"
#include "addToRunTimeSelectionTable.H"
#include "bandCompression.H"
#include "decompositionMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(boundaryFirstRenumber, 0);

    addToRunTimeSelectionTable
    (
        renumberMethod,
        boundaryFirstRenumber,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryFirstRenumber::boundaryFirstRenumber
(
    const dictionary& renumberDict
)
:
    renumberMethod(renumberDict),
    reverse_
    (
        renumberDict.found(typeName + "Coeffs")
      ? Switch(renumberDict.subDict(typeName + "Coeffs").lookup("reverse"))
      : Switch(true)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::boundaryFirstRenumber::renumber
(
    const polyMesh& mesh,
    const pointField& points
) const
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    // Distance to boundary. Minimum of neighbours.
    labelList distance(mesh.nCells(), -1);

    // New order
    labelList newOrder(mesh.nCells());
    label cellInOrder = 0;


    // Starting faces for walk. These are the zero distance faces.
    DynamicList<label> frontFaces(mesh.nFaces()-mesh.nInternalFaces());

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        //- Note: cannot check for empty since these are introduced by
        //  fvMeshSubset. Also most loops don't care about patch type.
        //if (!isA<emptyPolyPatch>(pp))
        {
            forAll(pp, i)
            {
                frontFaces.append(pp.start()+i);
            }
        }
    }

    // Loop over all frontFaces
    label currentDistance = 0;
    while (frontFaces.size() > 0)
    {
        DynamicList<label> frontCells(frontFaces.size());

        // Set all of frontFaces' neighbours to current distance

        forAll(frontFaces, i)
        {
            label faceI = frontFaces[i];

            if (mesh.isInternalFace(faceI))
            {
                label ownCellI = own[faceI];
                if (distance[ownCellI] == -1)
                {
                    distance[ownCellI] = currentDistance;
                    frontCells.append(ownCellI);
                }
                label neiCellI = nei[faceI];
                if (distance[neiCellI] == -1)
                {
                    distance[neiCellI] = currentDistance;
                    frontCells.append(neiCellI);
                }
            }
            else
            {
                label ownCellI = own[faceI];
                if (distance[ownCellI] == -1)
                {
                    distance[ownCellI] = currentDistance;
                    frontCells.append(ownCellI);
                }
            }
        }

        //Pout<< "For distance:" << currentDistance
        //    << " from " << frontFaces.size() << " faces to "
        //    << frontCells.size() << " cells." << endl;


        // TBD. Determine order within current shell (frontCells). For now
        // just add them.
        forAll(frontCells, i)
        {
            newOrder[cellInOrder] = frontCells[i];
            cellInOrder++;
        }

        // From cells to faces
        frontFaces.clear();
        forAll(frontCells, i)
        {
            label cellI = frontCells[i];
            const cell& cFaces = mesh.cells()[cellI];

            forAll(cFaces, i)
            {
                label faceI = cFaces[i];
                if (mesh.isInternalFace(faceI))
                {
                    label nbrCellI =
                    (
                        mesh.faceOwner()[faceI] == cellI
                      ? mesh.faceNeighbour()[faceI]
                      : mesh.faceOwner()[faceI]
                    );
                    if (distance[nbrCellI] == -1)
                    {
                        frontFaces.append(faceI);
                    }
                }
            }
        }

        currentDistance++;
    }

    // Return furthest away cell first
    if (reverse_)
    {
        reverse(newOrder);
    }

    //forAll(newOrder, i)
    //{
    //    label cellI = newOrder[i];
    //
    //    Pout<< "cell:" << cellI << endl;
    //    Pout<< " at distance:" << distance[cellI]
    //        << endl;
    //}

    return invert(newOrder.size(), newOrder);
}


// ************************************************************************* //
