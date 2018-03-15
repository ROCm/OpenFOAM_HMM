/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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

#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicFvMesh, 0);
    defineRunTimeSelectionTable(dynamicFvMesh, IOobject);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicFvMesh::dynamicFvMesh(const IOobject& io)
:
    fvMesh(io)
{}


Foam::dynamicFvMesh::dynamicFvMesh
(
    const IOobject& io,
    const zero,
    const bool syncPar
)
:
    fvMesh(io, Zero, syncPar)
{}


Foam::dynamicFvMesh::dynamicFvMesh
(
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    labelList&& allOwner,
    labelList&& allNeighbour,
    const bool syncPar
)
:
    fvMesh
    (
        io,
        std::move(points),
        std::move(faces),
        std::move(allOwner),
        std::move(allNeighbour),
        syncPar
    )
{}


Foam::dynamicFvMesh::dynamicFvMesh
(
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    cellList&& cells,
    const bool syncPar
)
:
    fvMesh
    (
        io,
        std::move(points),
        std::move(faces),
        std::move(cells),
        syncPar
    )
{}


// ************************************************************************* //
