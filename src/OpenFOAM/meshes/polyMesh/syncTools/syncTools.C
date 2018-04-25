/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

    label nBnd = mesh.nFaces()-mesh.nInternalFaces();

    neighbourCellData.setSize(nBnd);

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];
        const labelUList& faceCells = pp.faceCells();
        forAll(faceCells, i)
        {
            label bFacei = pp.start()+i-mesh.nInternalFaces();
            neighbourCellData[bFacei] = cellData[faceCells[i]];
        }
    }
    syncTools::swapBoundaryFacePositions(mesh, neighbourCellData);
}


Foam::bitSet Foam::syncTools::getMasterPoints(const polyMesh& mesh)
{
    bitSet isMasterPoint(mesh.nPoints());
    bitSet donePoint(mesh.nPoints());

    const globalMeshData& globalData = mesh.globalData();
    const labelList& meshPoints = globalData.coupledPatch().meshPoints();
    const labelListList& slaves = globalData.globalPointSlaves();
    const labelListList& transformedSlaves =
            globalData.globalPointTransformedSlaves();

    forAll(meshPoints, coupledPointi)
    {
        label meshPointi = meshPoints[coupledPointi];
        if
        (
            (
                slaves[coupledPointi].size()
              + transformedSlaves[coupledPointi].size()
            )
          > 0
        )
        {
            isMasterPoint.set(meshPointi);
        }
        donePoint.set(meshPointi);
    }


    // Do all other points
    // ~~~~~~~~~~~~~~~~~~~

    forAll(donePoint, pointi)
    {
        if (!donePoint.test(pointi))
        {
            isMasterPoint.set(pointi);
        }
    }

    return isMasterPoint;
}


Foam::bitSet Foam::syncTools::getMasterEdges(const polyMesh& mesh)
{
    bitSet isMasterEdge(mesh.nEdges());
    bitSet doneEdge(mesh.nEdges());

    const globalMeshData& globalData = mesh.globalData();
    const labelList& meshEdges = globalData.coupledPatchMeshEdges();
    const labelListList& slaves = globalData.globalEdgeSlaves();
    const labelListList& transformedSlaves =
        globalData.globalEdgeTransformedSlaves();

    forAll(meshEdges, coupledEdgeI)
    {
        label meshEdgeI = meshEdges[coupledEdgeI];
        if
        (
            (
                slaves[coupledEdgeI].size()
              + transformedSlaves[coupledEdgeI].size()
            )
          > 0
        )
        {
            isMasterEdge.set(meshEdgeI);
        }
        doneEdge.set(meshEdgeI);
    }


    // Do all other edges
    // ~~~~~~~~~~~~~~~~~~

    forAll(doneEdge, edgeI)
    {
        if (!doneEdge.test(edgeI))
        {
            isMasterEdge.set(edgeI);
        }
    }

    return isMasterEdge;
}


Foam::bitSet Foam::syncTools::getMasterFaces(const polyMesh& mesh)
{
    bitSet isMasterFace(mesh.nFaces(), true);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchi)
    {
        if (patches[patchi].coupled())
        {
            const coupledPolyPatch& pp =
                refCast<const coupledPolyPatch>(patches[patchi]);

            if (!pp.owner())
            {
                forAll(pp, i)
                {
                    isMasterFace.unset(pp.start()+i);
                }
            }
        }
    }

    return isMasterFace;
}


Foam::bitSet Foam::syncTools::getInternalOrMasterFaces
(
    const polyMesh& mesh
)
{
    bitSet isMasterFace(mesh.nFaces(), true);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled())
        {
            if (!refCast<const coupledPolyPatch>(pp).owner())
            {
                forAll(pp, i)
                {
                    isMasterFace.unset(pp.start()+i);
                }
            }
        }
        else
        {
            forAll(pp, i)
            {
                isMasterFace.unset(pp.start()+i);
            }
        }
    }

    return isMasterFace;
}


Foam::bitSet Foam::syncTools::getInternalOrCoupledFaces
(
    const polyMesh& mesh
)
{
    bitSet isMasterFace(mesh.nFaces(), true);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (!pp.coupled())
        {
            forAll(pp, i)
            {
                isMasterFace.unset(pp.start()+i);
            }
        }
    }

    return isMasterFace;
}


// ************************************************************************* //
