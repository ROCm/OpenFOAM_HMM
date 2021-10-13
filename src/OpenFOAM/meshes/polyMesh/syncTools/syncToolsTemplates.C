/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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
#include "polyMesh.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "globalMeshData.H"
#include "transformList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T, class CombineOp>
void Foam::syncTools::combine
(
    Map<T>& pointValues,
    const CombineOp& cop,
    const label index,
    const T& val
)
{
    auto iter = pointValues.find(index);

    if (iter.found())
    {
        cop(*iter, val);
    }
    else
    {
        pointValues.insert(index, val);
    }
}


template<class T, class CombineOp>
void Foam::syncTools::combine
(
    EdgeMap<T>& edgeValues,
    const CombineOp& cop,
    const edge& index,
    const T& val
)
{
    auto iter = edgeValues.find(index);

    if (iter.found())
    {
        cop(*iter, val);
    }
    else
    {
        edgeValues.insert(index, val);
    }
}


template<class T, class CombineOp, class TransformOp>
void Foam::syncTools::syncPointMap
(
    const polyMesh& mesh,
    Map<T>& pointValues,        // from mesh point label to value
    const CombineOp& cop,
    const TransformOp& top
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Synchronize multiple shared points.
    const globalMeshData& pd = mesh.globalData();

    // Values on shared points. Keyed on global shared index.
    Map<T> sharedPointValues(0);

    if (pd.nGlobalPoints() > 0)
    {
        // meshPoint per local index
        const labelList& sharedPtLabels = pd.sharedPointLabels();

        // global shared index per local index
        const labelList& sharedPtAddr = pd.sharedPointAddr();

        sharedPointValues.resize(sharedPtAddr.size());

        // Fill my entries in the shared points
        forAll(sharedPtLabels, i)
        {
            const auto fnd = pointValues.cfind(sharedPtLabels[i]);

            if (fnd.found())
            {
                combine
                (
                    sharedPointValues,
                    cop,
                    sharedPtAddr[i],    // index
                    *fnd                // value
                );
            }
        }
    }


    if (Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send
        for (const polyPatch& pp : patches)
        {
            const auto* ppp = isA<processorPolyPatch>(pp);

            if (ppp && pp.nPoints())
            {
                const auto& procPatch = *ppp;

                // Get data per patchPoint in neighbouring point numbers.

                const labelList& meshPts = procPatch.meshPoints();
                const labelList& nbrPts = procPatch.neighbPoints();

                // Extract local values. Create map from nbrPoint to value.
                // Note: how small initial size?
                Map<T> patchInfo(meshPts.size() / 20);

                forAll(meshPts, i)
                {
                    const auto iter = pointValues.cfind(meshPts[i]);

                    if (iter.found())
                    {
                        patchInfo.insert(nbrPts[i], *iter);
                    }
                }

                UOPstream toNeighb(procPatch.neighbProcNo(), pBufs);
                toNeighb << patchInfo;
            }
        }

        pBufs.finishedSends();

        // Receive and combine.
        for (const polyPatch& pp : patches)
        {
            const auto* ppp = isA<processorPolyPatch>(pp);

            if (ppp && pp.nPoints())
            {
                const auto& procPatch = *ppp;

                UIPstream fromNbr(procPatch.neighbProcNo(), pBufs);
                Map<T> nbrPatchInfo(fromNbr);

                // Transform
                top(procPatch, nbrPatchInfo);

                const labelList& meshPts = procPatch.meshPoints();

                // Only update those values which come from neighbour
                forAllConstIters(nbrPatchInfo, nbrIter)
                {
                    combine
                    (
                        pointValues,
                        cop,
                        meshPts[nbrIter.key()],
                        nbrIter.val()
                    );
                }
            }
        }
    }

    // Do the cyclics.
    for (const polyPatch& pp : patches)
    {
        const cyclicPolyPatch* cpp = isA<cyclicPolyPatch>(pp);

        if (cpp && cpp->owner())
        {
            // Owner does all.

            const cyclicPolyPatch& cycPatch = *cpp;
            const cyclicPolyPatch& nbrPatch = cycPatch.neighbPatch();

            const edgeList& coupledPoints = cycPatch.coupledPoints();
            const labelList& meshPtsA = cycPatch.meshPoints();
            const labelList& meshPtsB = nbrPatch.meshPoints();

            // Extract local values. Create map from coupled-edge to value.
            Map<T> half0Values(meshPtsA.size() / 20);
            Map<T> half1Values(half0Values.size());

            forAll(coupledPoints, i)
            {
                const edge& e = coupledPoints[i];

                const auto point0Fnd = pointValues.cfind(meshPtsA[e[0]]);

                if (point0Fnd.found())
                {
                    half0Values.insert(i, *point0Fnd);
                }

                const auto point1Fnd = pointValues.cfind(meshPtsB[e[1]]);

                if (point1Fnd.found())
                {
                    half1Values.insert(i, *point1Fnd);
                }
            }

            // Transform to receiving side
            top(cycPatch, half1Values);
            top(nbrPatch, half0Values);

            forAll(coupledPoints, i)
            {
                const edge& e = coupledPoints[i];

                const auto half0Fnd = half0Values.cfind(i);

                if (half0Fnd.found())
                {
                    combine
                    (
                        pointValues,
                        cop,
                        meshPtsB[e[1]],
                        *half0Fnd
                    );
                }

                const auto half1Fnd = half1Values.cfind(i);

                if (half1Fnd.found())
                {
                    combine
                    (
                        pointValues,
                        cop,
                        meshPtsA[e[0]],
                        *half1Fnd
                    );
                }
            }
        }
    }

    // Synchronize multiple shared points.
    if (pd.nGlobalPoints() > 0)
    {
        // meshPoint per local index
        const labelList& sharedPtLabels = pd.sharedPointLabels();
        // global shared index per local index
        const labelList& sharedPtAddr = pd.sharedPointAddr();

        // Reduce on master.

        if (Pstream::parRun())
        {
            if (Pstream::master())
            {
                // Receive the edges using shared points from the slave.
                for (const int slave : Pstream::subProcs())
                {
                    IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
                    Map<T> nbrValues(fromSlave);

                    // Merge neighbouring values with my values
                    forAllConstIters(nbrValues, iter)
                    {
                        combine
                        (
                            sharedPointValues,
                            cop,
                            iter.key(),     // edge
                            iter.val()      // value
                        );
                    }
                }

                // Send back
                for (const int slave : Pstream::subProcs())
                {
                    OPstream toSlave(Pstream::commsTypes::scheduled, slave);
                    toSlave << sharedPointValues;
                }
            }
            else
            {
                // Slave: send to master
                {
                    OPstream toMaster
                    (
                        Pstream::commsTypes::scheduled,
                        Pstream::masterNo()
                    );
                    toMaster << sharedPointValues;
                }
                // Receive merged values
                {
                    IPstream fromMaster
                    (
                        Pstream::commsTypes::scheduled,
                        Pstream::masterNo()
                    );
                    fromMaster >> sharedPointValues;
                }
            }
        }


        // Merge sharedPointValues (keyed on sharedPointAddr) into
        // pointValues (keyed on mesh points).

        // Map from global shared index to meshpoint
        Map<label> sharedToMeshPoint(2*sharedPtAddr.size());
        forAll(sharedPtAddr, i)
        {
            sharedToMeshPoint.insert(sharedPtAddr[i], sharedPtLabels[i]);
        }

        forAllConstIters(sharedToMeshPoint, iter)
        {
            // Do I have a value for my shared point
            const auto sharedFnd = sharedPointValues.cfind(iter.key());

            if (sharedFnd.found())
            {
                pointValues.set(*iter, *sharedFnd);
            }
        }
    }
}


template<class T, class CombineOp, class TransformOp>
void Foam::syncTools::syncEdgeMap
(
    const polyMesh& mesh,
    EdgeMap<T>& edgeValues,
    const CombineOp& cop,
    const TransformOp& top
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();


    // Do synchronisation without constructing globalEdge addressing
    // (since this constructs mesh edge addressing)


    // Swap proc patch info
    // ~~~~~~~~~~~~~~~~~~~~

    if (Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send
        for (const polyPatch& pp : patches)
        {
            const auto* ppp = isA<processorPolyPatch>(pp);

            if (ppp && pp.nEdges())
            {
                const auto& procPatch = *ppp;

                // Get data per patch edge in neighbouring edge.

                const edgeList& edges = procPatch.edges();
                const labelList& meshPts = procPatch.meshPoints();
                const labelList& nbrPts = procPatch.neighbPoints();

                EdgeMap<T> patchInfo(edges.size() / 20);

                for (const edge& e : edges)
                {
                    const edge meshEdge(meshPts[e[0]], meshPts[e[1]]);

                    const auto iter = edgeValues.cfind(meshEdge);

                    if (iter.found())
                    {
                        const edge nbrEdge(nbrPts[e[0]], nbrPts[e[1]]);
                        patchInfo.insert(nbrEdge, *iter);
                    }
                }

                UOPstream toNeighb(procPatch.neighbProcNo(), pBufs);
                toNeighb << patchInfo;
            }
        }

        pBufs.finishedSends();

        // Receive and combine.
        for (const polyPatch& pp : patches)
        {
            const auto* ppp = isA<processorPolyPatch>(pp);

            if (ppp && pp.nEdges())
            {
                const auto& procPatch = *ppp;

                EdgeMap<T> nbrPatchInfo;
                {
                    UIPstream fromNbr(procPatch.neighbProcNo(), pBufs);
                    fromNbr >> nbrPatchInfo;
                }

                // Apply transform to convert to this side properties.
                top(procPatch, nbrPatchInfo);


                // Only update those values which come from neighbour
                const labelList& meshPts = procPatch.meshPoints();

                forAllConstIters(nbrPatchInfo, nbrIter)
                {
                    const edge& e = nbrIter.key();
                    const edge meshEdge(meshPts[e[0]], meshPts[e[1]]);

                    combine
                    (
                        edgeValues,
                        cop,
                        meshEdge,           // edge
                        nbrIter.val()       // value
                    );
                }
            }
        }
    }


    // Swap cyclic info
    // ~~~~~~~~~~~~~~~~

    for (const polyPatch& pp : patches)
    {
        const cyclicPolyPatch* cpp = isA<cyclicPolyPatch>(pp);

        if (cpp && cpp->owner())
        {
            // Owner does all.

            const cyclicPolyPatch& cycPatch = *cpp;
            const cyclicPolyPatch& nbrPatch = cycPatch.neighbPatch();

            const edgeList& coupledEdges = cycPatch.coupledEdges();

            const labelList& meshPtsA = cycPatch.meshPoints();
            const edgeList& edgesA = cycPatch.edges();

            const labelList& meshPtsB = nbrPatch.meshPoints();
            const edgeList& edgesB = nbrPatch.edges();

            // Extract local values. Create map from edge to value.
            Map<T> half0Values(edgesA.size() / 20);
            Map<T> half1Values(half0Values.size());

            forAll(coupledEdges, edgei)
            {
                const edge& twoEdges = coupledEdges[edgei];

                {
                    const edge& e0 = edgesA[twoEdges[0]];
                    const edge meshEdge0(meshPtsA[e0[0]], meshPtsA[e0[1]]);

                    const auto iter = edgeValues.cfind(meshEdge0);

                    if (iter.found())
                    {
                        half0Values.insert(edgei, *iter);
                    }
                }
                {
                    const edge& e1 = edgesB[twoEdges[1]];
                    const edge meshEdge1(meshPtsB[e1[0]], meshPtsB[e1[1]]);

                    const auto iter = edgeValues.cfind(meshEdge1);

                    if (iter.found())
                    {
                        half1Values.insert(edgei, *iter);
                    }
                }
            }

            // Transform to this side
            top(cycPatch, half1Values);
            top(nbrPatch, half0Values);


            // Extract and combine information

            forAll(coupledEdges, edgei)
            {
                const edge& twoEdges = coupledEdges[edgei];

                const auto half1Fnd = half1Values.cfind(edgei);

                if (half1Fnd.found())
                {
                    const edge& e0 = edgesA[twoEdges[0]];
                    const edge meshEdge0(meshPtsA[e0[0]], meshPtsA[e0[1]]);

                    combine
                    (
                        edgeValues,
                        cop,
                        meshEdge0,      // edge
                        *half1Fnd       // value
                    );
                }

                const auto half0Fnd = half0Values.cfind(edgei);

                if (half0Fnd.found())
                {
                    const edge& e1 = edgesB[twoEdges[1]];
                    const edge meshEdge1(meshPtsB[e1[0]], meshPtsB[e1[1]]);

                    combine
                    (
                        edgeValues,
                        cop,
                        meshEdge1,      // edge
                        *half0Fnd       // value
                    );
                }
            }
        }
    }

    // Synchronize multiple shared points.
    // Problem is that we don't want to construct shared edges so basically
    // we do here like globalMeshData but then using sparse edge representation
    // (EdgeMap instead of mesh.edges())

    const globalMeshData& pd = mesh.globalData();
    const labelList& sharedPtAddr = pd.sharedPointAddr();
    const labelList& sharedPtLabels = pd.sharedPointLabels();

    // 1. Create map from meshPoint to globalShared index.
    Map<label> meshToShared(2*sharedPtLabels.size());
    forAll(sharedPtLabels, i)
    {
        meshToShared.insert(sharedPtLabels[i], sharedPtAddr[i]);
    }

    // Values on shared points. From two sharedPtAddr indices to a value.
    EdgeMap<T> sharedEdgeValues(meshToShared.size());

    // From shared edge to mesh edge. Used for merging later on.
    EdgeMap<edge> potentialSharedEdge(meshToShared.size());

    // 2. Find any edges using two global shared points. These will always be
    // on the outside of the mesh. (though might not be on coupled patch
    // if is single edge and not on coupled face)
    // Store value (if any) on sharedEdgeValues
    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); ++facei)
    {
        const face& f = mesh.faces()[facei];

        forAll(f, fp)
        {
            const label v0 = f[fp];
            const label v1 = f[f.fcIndex(fp)];

            const auto v0Fnd = meshToShared.cfind(v0);

            if (v0Fnd.found())
            {
                const auto v1Fnd = meshToShared.cfind(v1);

                if (v1Fnd.found())
                {
                    const edge meshEdge(v0, v1);

                    // edge in shared point labels
                    const edge sharedEdge(*v0Fnd, *v1Fnd);

                    // Store mesh edge as a potential shared edge.
                    potentialSharedEdge.insert(sharedEdge, meshEdge);

                    const auto edgeFnd = edgeValues.cfind(meshEdge);

                    if (edgeFnd.found())
                    {
                        // edge exists in edgeValues. See if already in map
                        // (since on same processor, e.g. cyclic)
                        combine
                        (
                            sharedEdgeValues,
                            cop,
                            sharedEdge, // edge
                            *edgeFnd    // value
                        );
                    }
                }
            }
        }
    }


    // Now sharedEdgeValues will contain per potential sharedEdge the value.
    // (potential since an edge having two shared points is not necessary a
    //  shared edge).
    // Reduce this on the master.

    if (Pstream::parRun())
    {
        if (Pstream::master())
        {
            // Receive the edges using shared points from the slave.
            for (const int slave : Pstream::subProcs())
            {
                IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
                EdgeMap<T> nbrValues(fromSlave);

                // Merge neighbouring values with my values
                forAllConstIters(nbrValues, iter)
                {
                    combine
                    (
                        sharedEdgeValues,
                        cop,
                        iter.key(),     // edge
                        iter.val()      // value
                    );
                }
            }

            // Send back
            for (const int slave : Pstream::subProcs())
            {
                OPstream toSlave(Pstream::commsTypes::scheduled, slave);
                toSlave << sharedEdgeValues;
            }
        }
        else
        {
            // Send to master
            {
                OPstream toMaster
                (
                    Pstream::commsTypes::scheduled,
                    Pstream::masterNo()
                );
                toMaster << sharedEdgeValues;
            }
            // Receive merged values
            {
                IPstream fromMaster
                (
                    Pstream::commsTypes::scheduled,
                    Pstream::masterNo()
                );
                fromMaster >> sharedEdgeValues;
            }
        }
    }


    // Merge sharedEdgeValues (keyed on sharedPointAddr) into edgeValues
    // (keyed on mesh points).

    // Loop over all my shared edges.
    forAllConstIters(potentialSharedEdge, iter)
    {
        const edge& sharedEdge = iter.key();
        const edge& meshEdge = iter.val();

        // Do I have a value for the shared edge?
        const auto sharedFnd = sharedEdgeValues.cfind(sharedEdge);

        if (sharedFnd.found())
        {
            combine
            (
                edgeValues,
                cop,
                meshEdge,       // edge
                *sharedFnd      // value
            );
        }
    }
}


template<class T, class CombineOp, class TransformOp>
void Foam::syncTools::syncPointList
(
    const polyMesh& mesh,
    List<T>& pointValues,
    const CombineOp& cop,
    const T& nullValue,
    const TransformOp& top
)
{
    if (pointValues.size() != mesh.nPoints())
    {
        FatalErrorInFunction
            << "Number of values " << pointValues.size()
            << " is not equal to the number of points in the mesh "
            << mesh.nPoints() << abort(FatalError);
    }

    mesh.globalData().syncPointData(pointValues, cop, top);
}


template<class T, class CombineOp, class TransformOp>
void Foam::syncTools::syncPointList
(
    const polyMesh& mesh,
    const labelUList& meshPoints,
    List<T>& pointValues,
    const CombineOp& cop,
    const T& nullValue,
    const TransformOp& top
)
{
    if (pointValues.size() != meshPoints.size())
    {
        FatalErrorInFunction
            << "Number of values " << pointValues.size()
            << " is not equal to the number of meshPoints "
            << meshPoints.size() << abort(FatalError);
    }
    const globalMeshData& gd = mesh.globalData();
    const indirectPrimitivePatch& cpp = gd.coupledPatch();
    const Map<label>& mpm = cpp.meshPointMap();

    List<T> cppFld(cpp.nPoints(), nullValue);

    forAll(meshPoints, i)
    {
        const auto iter = mpm.cfind(meshPoints[i]);

        if (iter.found())
        {
            cppFld[*iter] = pointValues[i];
        }
    }

    globalMeshData::syncData
    (
        cppFld,
        gd.globalPointSlaves(),
        gd.globalPointTransformedSlaves(),
        gd.globalPointSlavesMap(),
        gd.globalTransforms(),
        cop,
        top
    );

    forAll(meshPoints, i)
    {
        const auto iter = mpm.cfind(meshPoints[i]);

        if (iter.found())
        {
            pointValues[i] = cppFld[*iter];
        }
    }
}


template<class T, class CombineOp, class TransformOp, class FlipOp>
void Foam::syncTools::syncEdgeList
(
    const polyMesh& mesh,
    List<T>& edgeValues,
    const CombineOp& cop,
    const T& nullValue,
    const TransformOp& top,
    const FlipOp& fop
)
{
    if (edgeValues.size() != mesh.nEdges())
    {
        FatalErrorInFunction
            << "Number of values " << edgeValues.size()
            << " is not equal to the number of edges in the mesh "
            << mesh.nEdges() << abort(FatalError);
    }

    const edgeList& edges = mesh.edges();
    const globalMeshData& gd = mesh.globalData();
    const labelList& meshEdges = gd.coupledPatchMeshEdges();
    const indirectPrimitivePatch& cpp = gd.coupledPatch();
    const edgeList& cppEdges = cpp.edges();
    const labelList& mp = cpp.meshPoints();
    const globalIndexAndTransform& git = gd.globalTransforms();
    const mapDistribute& edgeMap = gd.globalEdgeSlavesMap();
    const bitSet& orientation = gd.globalEdgeOrientation();

    List<T> cppFld(meshEdges.size());
    forAll(meshEdges, i)
    {
        const edge& cppE = cppEdges[i];
        const label meshEdgei = meshEdges[i];
        const edge& meshE = edges[meshEdgei];

        // 1. is cpp edge oriented as mesh edge
        // 2. is cpp edge oriented same as master edge

        const int dir = edge::compare(meshE, edge(mp, cppE));
        if (dir == 0)
        {
            FatalErrorInFunction<< "Problem:"
                << " mesh edge:" << meshE.line(mesh.points())
                << " coupled edge:" << cppE.line(cpp.localPoints())
                << exit(FatalError);
        }

        const bool sameOrientation = ((dir == 1) == orientation[i]);

        if (sameOrientation)
        {
            cppFld[i] = edgeValues[meshEdgei];
        }
        else
        {
            cppFld[i] = fop(edgeValues[meshEdgei]);
        }
    }


    globalMeshData::syncData
    (
        cppFld,
        gd.globalEdgeSlaves(),
        gd.globalEdgeTransformedSlaves(),
        edgeMap,
        git,
        cop,
        top
    );

    // Extract back onto mesh
    forAll(meshEdges, i)
    {
        const edge& cppE = cppEdges[i];
        const label meshEdgei = meshEdges[i];
        const edge& meshE = edges[meshEdgei];

        // 1. is cpp edge oriented as mesh edge
        // 2. is cpp edge oriented same as master edge

        const int dir = edge::compare(meshE, edge(mp, cppE));
        const bool sameOrientation = ((dir == 1) == orientation[i]);

        if (sameOrientation)
        {
            edgeValues[meshEdges[i]] = cppFld[i];
        }
        else
        {
            edgeValues[meshEdges[i]] =  fop(cppFld[i]);
        }
    }
}


template<class T, class CombineOp, class TransformOp, class FlipOp>
void Foam::syncTools::syncEdgeList
(
    const polyMesh& mesh,
    const labelList& meshEdges,
    List<T>& edgeValues,
    const CombineOp& cop,
    const T& nullValue,
    const TransformOp& top,
    const FlipOp& fop
)
{
    if (edgeValues.size() != meshEdges.size())
    {
        FatalErrorInFunction
            << "Number of values " << edgeValues.size()
            << " is not equal to the number of meshEdges "
            << meshEdges.size() << abort(FatalError);
    }
    const edgeList& edges = mesh.edges();
    const globalMeshData& gd = mesh.globalData();
    const indirectPrimitivePatch& cpp = gd.coupledPatch();
    const edgeList& cppEdges = cpp.edges();
    const labelList& mp = cpp.meshPoints();
    const Map<label>& mpm = gd.coupledPatchMeshEdgeMap();
    const bitSet& orientation = gd.globalEdgeOrientation();

    List<T> cppFld(cpp.nEdges(), nullValue);

    forAll(meshEdges, i)
    {
        const label meshEdgei = meshEdges[i];
        const auto iter = mpm.cfind(meshEdgei);
        if (iter.found())
        {
            const label cppEdgei = iter();
            const edge& cppE = cppEdges[cppEdgei];
            const edge& meshE = edges[meshEdgei];

            // 1. is cpp edge oriented as mesh edge
            // 2. is cpp edge oriented same as master edge

            const int dir = edge::compare(meshE, edge(mp, cppE));
            if (dir == 0)
            {
                FatalErrorInFunction<< "Problem:"
                    << " mesh edge:" << meshE.line(mesh.points())
                    << " coupled edge:" << cppE.line(cpp.localPoints())
                    << exit(FatalError);
            }

            const bool sameOrientation = ((dir == 1) == orientation[i]);

            if (sameOrientation)
            {
                cppFld[cppEdgei] = edgeValues[i];
            }
            else
            {
                cppFld[cppEdgei] = fop(edgeValues[i]);
            }
        }
    }

    globalMeshData::syncData
    (
        cppFld,
        gd.globalEdgeSlaves(),
        gd.globalEdgeTransformedSlaves(),
        gd.globalEdgeSlavesMap(),
        gd.globalTransforms(),
        cop,
        top
    );

    forAll(meshEdges, i)
    {
        label meshEdgei = meshEdges[i];
        Map<label>::const_iterator iter = mpm.find(meshEdgei);
        if (iter != mpm.end())
        {
            label cppEdgei = iter();
            const edge& cppE = cppEdges[cppEdgei];
            const edge& meshE = edges[meshEdgei];

            const int dir = edge::compare(meshE, edge(mp, cppE));
            const bool sameOrientation = ((dir == 1) == orientation[i]);

            if (sameOrientation)
            {
                edgeValues[i] = cppFld[cppEdgei];
            }
            else
            {
                edgeValues[i] = fop(cppFld[cppEdgei]);
            }
        }
    }
}


template<class T, class CombineOp, class TransformOp>
void Foam::syncTools::syncBoundaryFaceList
(
    const polyMesh& mesh,
    UList<T>& faceValues,
    const CombineOp& cop,
    const TransformOp& top,
    const bool parRun
)
{
    // Offset (global to local) for start of boundaries
    const label boundaryOffset = mesh.nInternalFaces();

    if (faceValues.size() != mesh.nBoundaryFaces())
    {
        FatalErrorInFunction
            << "Number of values " << faceValues.size()
            << " is not equal to the number of boundary faces in the mesh "
            << mesh.nBoundaryFaces() << nl
            << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (parRun)
    {
        // Avoid mesh.globalData() - possible race condition

        if
        (
            is_contiguous<T>::value
         && Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
        )
        {
            const label startRequest = UPstream::nRequests();

            // Receive buffer
            List<T> receivedValues(mesh.nBoundaryFaces());

            // Set up reads
            for (const polyPatch& pp : patches)
            {
                const auto* ppp = isA<processorPolyPatch>(pp);

                if (ppp && pp.size())
                {
                    const auto& procPatch = *ppp;

                    SubList<T> fld
                    (
                        receivedValues,
                        pp.size(),
                        pp.start()-boundaryOffset
                    );

                    IPstream::read
                    (
                        Pstream::commsTypes::nonBlocking,
                        procPatch.neighbProcNo(),
                        fld.data_bytes(),
                        fld.size_bytes()
                    );
                }
            }

            // Set up writes
            for (const polyPatch& pp : patches)
            {
                const auto* ppp = isA<processorPolyPatch>(pp);

                if (ppp && pp.size())
                {
                    const auto& procPatch = *ppp;

                    const SubList<T> fld
                    (
                        faceValues,
                        pp.size(),
                        pp.start()-boundaryOffset
                    );

                    OPstream::write
                    (
                        Pstream::commsTypes::nonBlocking,
                        procPatch.neighbProcNo(),
                        fld.cdata_bytes(),
                        fld.size_bytes()
                    );
                }
            }

            // Wait for all comms to finish
            Pstream::waitRequests(startRequest);

            // Combine with existing data
            for (const polyPatch& pp : patches)
            {
                const auto* ppp = isA<processorPolyPatch>(pp);

                if (ppp && pp.size())
                {
                    const auto& procPatch = *ppp;

                    SubList<T> recvFld
                    (
                        receivedValues,
                        pp.size(),
                        pp.start()-boundaryOffset
                    );
                    const List<T>& fakeList = recvFld;
                    top(procPatch, const_cast<List<T>&>(fakeList));

                    SubList<T> patchValues
                    (
                        faceValues,
                        pp.size(),
                        pp.start()-boundaryOffset
                    );

                    forAll(patchValues, i)
                    {
                        cop(patchValues[i], recvFld[i]);
                    }
                }
            }
        }
        else
        {
            PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

            // Send
            for (const polyPatch& pp : patches)
            {
                const auto* ppp = isA<processorPolyPatch>(pp);

                if (ppp && pp.size())
                {
                    const auto& procPatch = *ppp;

                    const SubList<T> fld
                    (
                        faceValues,
                        pp.size(),
                        pp.start()-boundaryOffset
                    );

                    UOPstream toNbr(procPatch.neighbProcNo(), pBufs);
                    toNbr << fld;;
                }
            }

            pBufs.finishedSends();

            // Receive and combine.
            for (const polyPatch& pp : patches)
            {
                const auto* ppp = isA<processorPolyPatch>(pp);

                if (ppp && pp.size())
                {
                    const auto& procPatch = *ppp;

                    List<T> recvFld(pp.size());

                    UIPstream fromNbr(procPatch.neighbProcNo(), pBufs);
                    fromNbr >> recvFld;

                    top(procPatch, recvFld);

                    SubList<T> patchValues
                    (
                        faceValues,
                        pp.size(),
                        pp.start()-boundaryOffset
                    );

                    forAll(patchValues, i)
                    {
                        cop(patchValues[i], recvFld[i]);
                    }
                }
            }
        }
    }

    // Do the cyclics.
    for (const polyPatch& pp : patches)
    {
        const cyclicPolyPatch* cpp = isA<cyclicPolyPatch>(pp);

        if (cpp && cpp->owner())
        {
            // Owner does all.

            const cyclicPolyPatch& cycPatch = *cpp;
            const cyclicPolyPatch& nbrPatch = cycPatch.neighbPatch();
            const label patchSize = cycPatch.size();

            SubList<T> ownPatchValues
            (
                faceValues,
                patchSize,
                cycPatch.start()-boundaryOffset
            );

            SubList<T> nbrPatchValues
            (
                faceValues,
                patchSize,
                nbrPatch.start()-boundaryOffset
            );

            // Transform (copy of) data on both sides
            List<T> ownVals(ownPatchValues);
            top(nbrPatch, ownVals);

            List<T> nbrVals(nbrPatchValues);
            top(cycPatch, nbrVals);

            forAll(ownPatchValues, i)
            {
                cop(ownPatchValues[i], nbrVals[i]);
            }

            forAll(nbrPatchValues, i)
            {
                cop(nbrPatchValues[i], ownVals[i]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<unsigned Width, class CombineOp>
void Foam::syncTools::syncFaceList
(
    const polyMesh& mesh,
    const bool isBoundaryOnly,
    PackedList<Width>& faceValues,
    const CombineOp& cop,
    const bool parRun
)
{
    // Offset (global to local) for start of boundaries
    const label boundaryOffset = (isBoundaryOnly ? mesh.nInternalFaces() : 0);

    // Check size
    if (faceValues.size() != (mesh.nFaces() - boundaryOffset))
    {
        FatalErrorInFunction
            << "Number of values " << faceValues.size()
            << " is not equal to the number of "
            << (isBoundaryOnly ? "boundary" : "mesh") << " faces "
            << ((mesh.nFaces() - boundaryOffset)) << nl
            << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (parRun)
    {
        const label startRequest = UPstream::nRequests();

        // Receive buffers
        PtrList<PackedList<Width>> recvInfos(patches.size());

        // Set up reads
        for (const polyPatch& pp : patches)
        {
            const auto* ppp = isA<processorPolyPatch>(pp);

            if (ppp && pp.size())
            {
                const auto& procPatch = *ppp;
                const label patchi = pp.index();
                const label patchSize = pp.size();

                recvInfos.set(patchi, new PackedList<Width>(patchSize));
                PackedList<Width>& recvInfo = recvInfos[patchi];

                IPstream::read
                (
                    Pstream::commsTypes::nonBlocking,
                    procPatch.neighbProcNo(),
                    recvInfo.data_bytes(),
                    recvInfo.size_bytes()
                );
            }
        }

        // Send buffers
        PtrList<PackedList<Width>> sendInfos(patches.size());

        // Set up writes
        for (const polyPatch& pp : patches)
        {
            const auto* ppp = isA<processorPolyPatch>(pp);

            if (ppp && pp.size())
            {
                const auto& procPatch = *ppp;
                const label patchi = pp.index();

                const labelRange range
                (
                    pp.start()-boundaryOffset,
                    pp.size()
                );
                sendInfos.set
                (
                    patchi,
                    new PackedList<Width>(faceValues, range)
                );
                PackedList<Width>& sendInfo = sendInfos[patchi];

                OPstream::write
                (
                    Pstream::commsTypes::nonBlocking,
                    procPatch.neighbProcNo(),
                    sendInfo.cdata_bytes(),
                    sendInfo.size_bytes()
                );
            }
        }

        // Wait for all comms to finish
        Pstream::waitRequests(startRequest);

        // Combine with existing data
        for (const polyPatch& pp : patches)
        {
            const auto* ppp = isA<processorPolyPatch>(pp);

            if (ppp && pp.size())
            {
                const label patchi = pp.index();
                const label patchSize = pp.size();

                const PackedList<Width>& recvInfo = recvInfos[patchi];

                // Combine (bitwise)
                label bFacei = pp.start()-boundaryOffset;
                for (label i = 0; i < patchSize; ++i)
                {
                    unsigned int recvVal = recvInfo[i];
                    unsigned int faceVal = faceValues[bFacei];

                    cop(faceVal, recvVal);
                    faceValues.set(bFacei, faceVal);

                    ++bFacei;
                }
            }
        }
    }


    // Do the cyclics.
    for (const polyPatch& pp : patches)
    {
        const cyclicPolyPatch* cpp = isA<cyclicPolyPatch>(pp);

        if (cpp && cpp->owner())
        {
            // Owner does all.

            const cyclicPolyPatch& cycPatch = *cpp;
            const cyclicPolyPatch& nbrPatch = cycPatch.neighbPatch();
            const label patchSize = cycPatch.size();

            label face0 = cycPatch.start()-boundaryOffset;
            label face1 = nbrPatch.start()-boundaryOffset;
            for (label i = 0; i < patchSize; ++i)
            {
                unsigned int val0 = faceValues[face0];
                unsigned int val1 = faceValues[face1];

                unsigned int t = val0;
                cop(t, val1);
                faceValues[face0] = t;

                cop(val1, val0);
                faceValues[face1] = val1;

                ++face0;
                ++face1;
            }
        }
    }
}


template<class T>
void Foam::syncTools::swapBoundaryCellList
(
    const polyMesh& mesh,
    const UList<T>& cellData,
    List<T>& neighbourCellData
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
        label bFacei = pp.offset();

        for (const label celli : pp.faceCells())
        {
            neighbourCellData[bFacei] = cellData[celli];
            ++bFacei;
        }
    }
    syncTools::swapBoundaryFaceList(mesh, neighbourCellData);
}


template<unsigned Width, class CombineOp>
void Foam::syncTools::syncFaceList
(
    const polyMesh& mesh,
    PackedList<Width>& faceValues,
    const CombineOp& cop,
    const bool parRun
)
{
    syncFaceList(mesh, false, faceValues, cop, parRun);
}


template<unsigned Width, class CombineOp>
void Foam::syncTools::syncBoundaryFaceList
(
    const polyMesh& mesh,
    PackedList<Width>& faceValues,
    const CombineOp& cop,
    const bool parRun
)
{
    syncFaceList(mesh, true, faceValues, cop, parRun);
}


template<unsigned Width>
void Foam::syncTools::swapFaceList
(
    const polyMesh& mesh,
    PackedList<Width>& faceValues
)
{
    syncFaceList(mesh, faceValues, eqOp<unsigned int>());
}


template<unsigned Width>
void Foam::syncTools::swapBoundaryFaceList
(
    const polyMesh& mesh,
    PackedList<Width>& faceValues
)
{
    syncBoundaryFaceList(mesh, faceValues, eqOp<unsigned int>());
}


template<unsigned Width, class CombineOp>
void Foam::syncTools::syncPointList
(
    const polyMesh& mesh,
    PackedList<Width>& pointValues,
    const CombineOp& cop,
    const unsigned int nullValue
)
{
    if (pointValues.size() != mesh.nPoints())
    {
        FatalErrorInFunction
            << "Number of values " << pointValues.size()
            << " is not equal to the number of points in the mesh "
            << mesh.nPoints() << abort(FatalError);
    }

    const globalMeshData& gd = mesh.globalData();
    const labelList& meshPoints = gd.coupledPatch().meshPoints();

    List<unsigned int> cppFld(gd.globalPointSlavesMap().constructSize());
    forAll(meshPoints, i)
    {
        cppFld[i] = pointValues[meshPoints[i]];
    }

    globalMeshData::syncData
    (
        cppFld,
        gd.globalPointSlaves(),
        gd.globalPointTransformedSlaves(),
        gd.globalPointSlavesMap(),
        cop
    );

    // Extract back to mesh
    forAll(meshPoints, i)
    {
        pointValues[meshPoints[i]] = cppFld[i];
    }
}


template<unsigned Width, class CombineOp>
void Foam::syncTools::syncEdgeList
(
    const polyMesh& mesh,
    PackedList<Width>& edgeValues,
    const CombineOp& cop,
    const unsigned int nullValue
)
{
    if (edgeValues.size() != mesh.nEdges())
    {
        FatalErrorInFunction
            << "Number of values " << edgeValues.size()
            << " is not equal to the number of edges in the mesh "
            << mesh.nEdges() << abort(FatalError);
    }

    const globalMeshData& gd = mesh.globalData();
    const labelList& meshEdges = gd.coupledPatchMeshEdges();

    List<unsigned int> cppFld(gd.globalEdgeSlavesMap().constructSize());
    forAll(meshEdges, i)
    {
        cppFld[i] = edgeValues[meshEdges[i]];
    }

    globalMeshData::syncData
    (
        cppFld,
        gd.globalEdgeSlaves(),
        gd.globalEdgeTransformedSlaves(),
        gd.globalEdgeSlavesMap(),
        cop
    );

    // Extract back to mesh
    forAll(meshEdges, i)
    {
        edgeValues[meshEdges[i]] = cppFld[i];
    }
}


// ************************************************************************* //
