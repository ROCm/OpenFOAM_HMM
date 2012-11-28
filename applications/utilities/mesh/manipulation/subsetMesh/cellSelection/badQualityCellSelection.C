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

#include "badQualityCellSelection.H"
#include "addToRunTimeSelectionTable.H"
#include "faceSet.H"
#include "polyMesh.H"
#include "motionSmoother.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cellSelections
{
    defineTypeNameAndDebug(badQualityCellSelection, 0);
    addToRunTimeSelectionTable
    (
        cellSelection,
        badQualityCellSelection,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellSelections::badQualityCellSelection::badQualityCellSelection
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    cellSelection(name, mesh, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellSelections::badQualityCellSelection::~badQualityCellSelection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cellSelections::badQualityCellSelection::select
(
    boolList& selectedCell
) const
{
    //- Delete cell of any face in error
    faceSet faces(mesh_, "meshQualityFaces", mesh_.nFaces()/100+1);
    motionSmoother::checkMesh(false, mesh_, dict_, faces);
    label nFaces = returnReduce(faces.size(), sumOp<label>());
    if (nFaces > 0)
    {
        faces.sync(mesh_);
        forAllConstIter(faceSet, faces, iter)
        {
            label faceI = iter.key();
            selectedCell[mesh_.faceOwner()[faceI]] = false;
            if (mesh_.isInternalFace(faceI))
            {
                selectedCell[mesh_.faceNeighbour()[faceI]] = false;
            }
        }
    }
}


// ************************************************************************* //
