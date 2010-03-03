/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

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


    // Do multiple shared points. Min. proc is master
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const labelList& sharedPointAddr =
        mesh.globalData().sharedPointAddr();

    labelList minProc(mesh.globalData().nGlobalPoints(), labelMax);

    UIndirectList<label>(minProc, sharedPointAddr) = Pstream::myProcNo();

    Pstream::listCombineGather(minProc, minEqOp<label>());
    Pstream::listCombineScatter(minProc);

    const labelList& sharedPointLabels =
        mesh.globalData().sharedPointLabels();

    forAll(sharedPointAddr, i)
    {
        if (minProc[sharedPointAddr[i]] == Pstream::myProcNo())
        {
            isMasterPoint.set(sharedPointLabels[i], 1u);
        }
        donePoint.set(sharedPointLabels[i], 1u);
    }


    // Do other points on coupled patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchI)
    {
        if (patches[patchI].coupled())
        {
            const coupledPolyPatch& pp =
                refCast<const coupledPolyPatch>(patches[patchI]);

            const labelList& meshPoints = pp.meshPoints();

            forAll(meshPoints, i)
            {
                label pointI = meshPoints[i];

                if (donePoint.get(pointI) == 0u)
                {
                    donePoint.set(pointI, 1u);

                    if (pp.owner())
                    {
                        isMasterPoint.set(pointI, 1u);
                    }
                }
            }
        }
    }


    // Do all other points
    // ~~~~~~~~~~~~~~~~~~~

    forAll(donePoint, pointI)
    {
        if (donePoint.get(pointI) == 0u)
        {
            donePoint.set(pointI, 1u);
            isMasterPoint.set(pointI, 1u);
        }
    }

    return isMasterPoint;
}


// Determines for every edge whether it is coupled and if so sets only one.
Foam::PackedBoolList Foam::syncTools::getMasterEdges(const polyMesh& mesh)
{
    PackedBoolList isMasterEdge(mesh.nEdges());
    PackedBoolList doneEdge(mesh.nEdges());


    // Do multiple shared edges. Min. proc is master
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const labelList& sharedEdgeAddr =
        mesh.globalData().sharedEdgeAddr();

    labelList minProc(mesh.globalData().nGlobalEdges(), labelMax);

    UIndirectList<label>(minProc, sharedEdgeAddr) = Pstream::myProcNo();

    Pstream::listCombineGather(minProc, minEqOp<label>());
    Pstream::listCombineScatter(minProc);

    const labelList& sharedEdgeLabels =
        mesh.globalData().sharedEdgeLabels();

    forAll(sharedEdgeAddr, i)
    {
        if (minProc[sharedEdgeAddr[i]] == Pstream::myProcNo())
        {
            isMasterEdge.set(sharedEdgeLabels[i], 1u);
        }
        doneEdge.set(sharedEdgeLabels[i], 1u);
    }


    // Do other edges on coupled patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchI)
    {
        if (patches[patchI].coupled())
        {
            const coupledPolyPatch& pp =
                refCast<const coupledPolyPatch>(patches[patchI]);

            const labelList& meshEdges = pp.meshEdges();

            forAll(meshEdges, i)
            {
                label edgeI = meshEdges[i];

                if (doneEdge.get(edgeI) == 0u)
                {
                    doneEdge.set(edgeI, 1u);

                    if (pp.owner())
                    {
                        isMasterEdge.set(edgeI, 1u);
                    }
                }
            }
        }
    }


    // Do all other edges
    // ~~~~~~~~~~~~~~~~~~

    forAll(doneEdge, edgeI)
    {
        if (doneEdge.get(edgeI) == 0u)
        {
            doneEdge.set(edgeI, 1u);
            isMasterEdge.set(edgeI, 1u);
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
