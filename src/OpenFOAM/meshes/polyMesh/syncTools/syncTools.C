/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "syncTools.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
void Foam::syncTools::transform::operator()
(
    const coupledPolyPatch&,
    Field<label>&
) const
{}
template<>
void Foam::syncTools::transform::operator()
(
    const coupledPolyPatch&,
    Map<label>&
) const
{}
template<>
void Foam::syncTools::transform::operator()
(
    const coupledPolyPatch&,
    EdgeMap<label>&
) const
{}


template<>
void Foam::syncTools::transform::operator()
(
    const coupledPolyPatch&,
    Field<scalar>&
) const
{}
template<>
void Foam::syncTools::transform::operator()
(
    const coupledPolyPatch&,
    Map<scalar>&
) const
{}
template<>
void Foam::syncTools::transform::operator()
(
    const coupledPolyPatch&,
    EdgeMap<scalar>&
) const
{}


template<>
void Foam::syncTools::transform::operator()
(
    const coupledPolyPatch&,
    Field<bool>&
) const
{}
template<>
void Foam::syncTools::transform::operator()
(
    const coupledPolyPatch&,
    Map<bool>&
) const
{}
template<>
void Foam::syncTools::transform::operator()
(
    const coupledPolyPatch&,
    EdgeMap<bool>&
) const
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Does anyone have couples? Since meshes might have 0 cells and 0 proc
// boundaries need to reduce this info.
bool Foam::syncTools::hasCouples(const polyBoundaryMesh& patches)
{
    bool hasAnyCouples = false;

    forAll(patches, patchI)
    {
        if (patches[patchI].coupled())
        {
            hasAnyCouples = true;
            break;
        }
    }
    return returnReduce(hasAnyCouples, orOp<bool>());
}


// Determines for every point whether it is coupled and if so sets only one.
Foam::PackedBoolList Foam::syncTools::getMasterPoints(const polyMesh& mesh)
{
    PackedBoolList isMasterPoint(mesh.nPoints());
    PackedBoolList donePoint(mesh.nPoints());

    const globalMeshData& globalData = mesh.globalData();
    const labelList& meshPoints = globalData.coupledPatch().meshPoints();
    const labelListList& pointSlaves = globalData.globalPointAllSlaves();

    forAll(meshPoints, coupledPointI)
    {
        label meshPointI = meshPoints[coupledPointI];
        if (pointSlaves[coupledPointI].size() > 0)
        {
            isMasterPoint[meshPointI] = true;
        }
        donePoint[meshPointI] = true;
    }


    // Do all other points
    // ~~~~~~~~~~~~~~~~~~~

    forAll(donePoint, pointI)
    {
        if (!donePoint[pointI])
        {
            isMasterPoint[pointI] = true;
        }
    }

    return isMasterPoint;
}


// Determines for every edge whether it is coupled and if so sets only one.
Foam::PackedBoolList Foam::syncTools::getMasterEdges(const polyMesh& mesh)
{
    PackedBoolList isMasterEdge(mesh.nEdges());
    PackedBoolList doneEdge(mesh.nEdges());

    const globalMeshData& globalData = mesh.globalData();
    const labelList& meshEdges = globalData.coupledPatchMeshEdges();
    const labelListList& edgeSlaves = globalData.globalEdgeAllSlaves();

    forAll(meshEdges, coupledEdgeI)
    {
        label meshEdgeI = meshEdges[coupledEdgeI];
        if (edgeSlaves[coupledEdgeI].size() > 0)
        {
            isMasterEdge[meshEdgeI] = true;
        }
        doneEdge[meshEdgeI] = true;
    }


    // Do all other edges
    // ~~~~~~~~~~~~~~~~~~

    forAll(doneEdge, edgeI)
    {
        if (!doneEdge[edgeI])
        {
            isMasterEdge[edgeI] = true;
        }
    }

    return isMasterEdge;
}


// Determines for every face whether it is coupled and if so sets only one.
Foam::PackedBoolList Foam::syncTools::getMasterFaces(const polyMesh& mesh)
{
    PackedBoolList isMasterFace(mesh.nFaces(), 1);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchI)
    {
        if (patches[patchI].coupled())
        {
            const coupledPolyPatch& pp =
                refCast<const coupledPolyPatch>(patches[patchI]);

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


// ************************************************************************* //
