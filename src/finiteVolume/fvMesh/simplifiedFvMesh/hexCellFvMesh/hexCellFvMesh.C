/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "hexCellFvMesh.H"
#include "emptyPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

namespace Foam
{
namespace simplifiedMeshes
{
    defineTypeNameAndDebug(hexCellFvMesh, 0);

    addToRunTimeSelectionTable
    (
        simplifiedFvMesh,
        hexCellFvMesh,
        time
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simplifiedMeshes::hexCellFvMesh::hexCellFvMesh
(
    const Time& runTime,
    const scalar d
)
:
    simplifiedFvMesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        pointField
        (
            {
                point(0, 0, 0),
                point(d, 0, 0),
                point(d, d, 0),
                point(0, d, 0),
                point(0, 0, d),
                point(d, 0, d),
                point(d, d, d),
                point(0, d, d)
            }
        ),
        faceList(cellModel::ref(cellModel::HEX).modelFaces()),
        labelList(6, Zero),
        labelList()
    )
{
    List<polyPatch*> patches(1);

    patches[0] = new emptyPolyPatch
    (
        "boundary",
        6,
        0,
        0,
        boundaryMesh(),
        emptyPolyPatch::typeName
    );

    addFvPatches(patches);

    // Note: no sets or zones created
}


// ************************************************************************* //
