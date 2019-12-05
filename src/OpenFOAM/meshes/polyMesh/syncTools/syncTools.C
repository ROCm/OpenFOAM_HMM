/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "syncTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::syncTools::swapBoundaryCellPositions
(
    const polyMesh& mesh,
    const UList<point>& cellData,
    List<point>& neighbourCellData
)
{
    if (cellData.size() != mesh.nCells())
    {
        FatalErrorInFunction
            << "Number of cell values " << cellData.size()
            << " is not equal to the number of cells in the mesh "
            << mesh.nCells() << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    neighbourCellData.resize(mesh.nBoundaryFaces());

    for (const polyPatch& pp : patches)
    {
        label bFacei = pp.start()-mesh.nInternalFaces();

        const labelUList& faceCells = pp.faceCells();

        for (const label celli : faceCells)
        {
            neighbourCellData[bFacei] = cellData[celli];
            ++bFacei;
        }
    }
    syncTools::swapBoundaryFacePositions(mesh, neighbourCellData);
}


Foam::bitSet Foam::syncTools::getMasterPoints(const polyMesh& mesh)
{
    bitSet isMaster(mesh.nPoints());
    bitSet unvisited(mesh.nPoints(), true);

    const globalMeshData& globalData = mesh.globalData();
    const labelList& meshPoints = globalData.coupledPatch().meshPoints();
    const labelListList& slaves = globalData.globalPointSlaves();
    const labelListList& transformedSlaves =
            globalData.globalPointTransformedSlaves();

    forAll(meshPoints, i)
    {
        const label meshPointi = meshPoints[i];

        if (!slaves[i].empty() || !transformedSlaves[i].empty())
        {
            isMaster.set(meshPointi);
        }
        unvisited.unset(meshPointi);
    }

    // Add in all unvisited points
    isMaster |= unvisited;

    return isMaster;
}


Foam::bitSet Foam::syncTools::getMasterEdges(const polyMesh& mesh)
{
    bitSet isMaster(mesh.nEdges());
    bitSet unvisited(mesh.nEdges(), true);

    const globalMeshData& globalData = mesh.globalData();
    const labelList& meshEdges = globalData.coupledPatchMeshEdges();
    const labelListList& slaves = globalData.globalEdgeSlaves();
    const labelListList& transformedSlaves =
        globalData.globalEdgeTransformedSlaves();

    forAll(meshEdges, i)
    {
        const label meshEdgei = meshEdges[i];

        if (!slaves[i].empty() || !transformedSlaves[i].empty())
        {
            isMaster.set(meshEdgei);
        }
        unvisited.unset(meshEdgei);
    }

    // Add in all unvisited edges
    isMaster |= unvisited;

    return isMaster;
}


Foam::bitSet Foam::syncTools::getMasterFaces(const polyMesh& mesh)
{
    bitSet isMaster(mesh.nFaces(), true);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    for (const polyPatch& pp : patches)
    {
        if (pp.coupled())
        {
            if (!refCast<const coupledPolyPatch>(pp).owner())
            {
                isMaster.unset(pp.range());
            }
        }
    }

    return isMaster;
}


Foam::bitSet Foam::syncTools::getInternalOrMasterFaces
(
    const polyMesh& mesh
)
{
    bitSet isMaster(mesh.nFaces(), true);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    for (const polyPatch& pp : patches)
    {
        if (pp.coupled())
        {
            if (!refCast<const coupledPolyPatch>(pp).owner())
            {
                isMaster.unset(pp.range());
            }
        }
        else
        {
            isMaster.unset(pp.range());
        }
    }

    return isMaster;
}


Foam::bitSet Foam::syncTools::getInternalOrCoupledFaces
(
    const polyMesh& mesh
)
{
    bitSet isMaster(mesh.nFaces(), true);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    for (const polyPatch& pp : patches)
    {
        if (!pp.coupled())
        {
            isMaster.unset(pp.range());
        }
    }

    return isMaster;
}


// ************************************************************************* //
