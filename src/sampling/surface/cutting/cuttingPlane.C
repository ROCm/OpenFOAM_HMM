/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

#include "cuttingPlane.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::cuttingPlane::debug(Foam::debug::debugSwitch("cuttingPlane", 0));


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cuttingPlane::cuttingPlane(const plane& pln)
:
    plane(pln)
{}


Foam::cuttingPlane::cuttingPlane
(
    const plane& pln,
    const primitiveMesh& mesh,
    const bool triangulate,
    const bitSet& cellIdLabels
)
:
    plane(pln)
{
    performCut(mesh, triangulate, cellIdLabels);
}


Foam::cuttingPlane::cuttingPlane
(
    const plane& pln,
    const primitiveMesh& mesh,
    const bool triangulate,
    bitSet&& cellIdLabels
)
:
    plane(pln)
{
    performCut(mesh, triangulate, cellIdLabels);
}


Foam::cuttingPlane::cuttingPlane
(
    const plane& pln,
    const primitiveMesh& mesh,
    const bool triangulate,
    const labelUList& cellIdLabels
)
:
    plane(pln)
{
    performCut(mesh, triangulate, cellIdLabels);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cuttingPlane::performCut
(
    const primitiveMesh& mesh,
    const bool triangulate,
    bitSet&& cellIdLabels
)
{
    MeshStorage::clear();
    meshCells_.clear();

    // Pre-populate with restriction
    bitSet cellCuts(std::move(cellIdLabels));

    if (cellCuts.size())
    {
        cellCuts.resize(mesh.nCells());
    }

    // For each mesh point, the encoded side (0,1,2) of the plane
    PackedList<2> sides;

    // Determine cells that are (likely) cut
    // - some ambiguity when plane is exactly between cells
    const label nFaceCuts = calcCellCuts(mesh, sides, cellCuts);

    // Find closed loop from cell cuts
    walkCellCuts(mesh, cellCuts, sides, triangulate, nFaceCuts);
}


void Foam::cuttingPlane::performCut
(
    const primitiveMesh& mesh,
    const bool triangulate,
    const bitSet& cellIdLabels
)
{
    bitSet subsetCells(cellIdLabels);

    performCut(mesh, triangulate, std::move(subsetCells));
}


void Foam::cuttingPlane::performCut
(
    const primitiveMesh& mesh,
    const bool triangulate,
    const labelUList& cellIdLabels
)
{
    bitSet subsetCells;

    if (notNull(cellIdLabels))
    {
        // Pre-populate with restriction
        subsetCells.resize(mesh.nCells());
        subsetCells.set(cellIdLabels);
    }

    performCut(mesh, triangulate, std::move(subsetCells));
}


void Foam::cuttingPlane::remapFaces(const labelUList& faceMap)
{
    if (notNull(faceMap) && !faceMap.empty())
    {
        MeshStorage::remapFaces(faceMap);

        List<label> remappedCells(faceMap.size());
        forAll(faceMap, facei)
        {
            remappedCells[facei] = meshCells_[faceMap[facei]];
        }
        meshCells_.transfer(remappedCells);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::cuttingPlane::operator=(const cuttingPlane& rhs)
{
    if (this == &rhs)
    {
        FatalErrorInFunction
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    static_cast<MeshStorage&>(*this) = rhs;
    static_cast<plane&>(*this) = rhs;
    meshCells_ = rhs.meshCells();
}


// ************************************************************************* //
