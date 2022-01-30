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

#include "cellBitSet.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellBitSet, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellBitSet::cellBitSet(const polyMesh& mesh)
:
    cellBitSet(mesh, false)
{}


Foam::cellBitSet::cellBitSet(const polyMesh& mesh, const bool val)
:
    topoBitSet(mesh, "cellBitSet", mesh.nCells(), val)
{}


Foam::cellBitSet::cellBitSet
(
    const polyMesh& mesh,
    const bitSet& bits
)
:
    topoBitSet(mesh, "cellBitSet", mesh.nCells(), bits)
{}


Foam::cellBitSet::cellBitSet
(
    const polyMesh& mesh,
    bitSet&& bits
)
:
    topoBitSet(mesh, "cellBitSet", mesh.nCells(), std::move(bits))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cellBitSet::maxSize(const polyMesh& mesh) const
{
    return mesh.nCells();
}


void Foam::cellBitSet::writeDebug
(
    Ostream& os,
    const primitiveMesh& mesh,
    const label maxLen
) const
{
    topoSet::writeDebug(os, mesh.cellCentres(), maxLen);
}


// ************************************************************************* //
