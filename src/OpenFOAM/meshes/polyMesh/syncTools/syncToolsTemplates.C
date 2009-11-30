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
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "globalMeshData.H"
#include "contiguous.H"
#include "transform.H"
#include "transformList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Combine val with existing value at index
template <class T, class CombineOp>
void Foam::syncTools::combine
(
    Map<T>& pointValues,
    const CombineOp& cop,
    const label index,
    const T& val
)
{
    typename Map<T>::iterator iter = pointValues.find(index);

    if (iter != pointValues.end())
    {
        cop(iter(), val);
    }
    else
    {
        pointValues.insert(index, val);
    }
}


// Combine val with existing value at (implicit index) e.
template <class T, class CombineOp>
void Foam::syncTools::combine
(
    EdgeMap<T>& edgeValues,
    const CombineOp& cop,
    const edge& index,
    const T& val
)
{
    typename EdgeMap<T>::iterator iter = edgeValues.find(index);

    if (iter != edgeValues.end())
    {
        cop(iter(), val);
    }
    else
    {
        edgeValues.insert(index, val);
    }
}


template <class T, class CombineOp, class TransformOp>
void Foam::syncTools::syncPointMap
(
    const polyMesh& mesh,
    Map<T>& pointValues,        // from mesh point label to value
    const CombineOp& cop,
    const TransformOp& top
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (!hasCouples(patches))
    {
        return;
    }

    if (Pstream::parRun())
    {
        // Send

        forAll(patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].nPoints() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                //TBD: optimisation for single subpatch.

                List<Map<T> > patchInfo(procPatch.patchIDs().size());

                forAll(procPatch.patchIDs(), subI)
                {
                    label subPatchI = procPatch.patchIDs()[subI];

                    const labelList& subMeshPts =
                        procPatch.subMeshPoints()[subI];

                    Map<T>& subPatchInfo = patchInfo[subI];
                    subPatchInfo.resize(subMeshPts.size()/20);
                    forAll(subMeshPts, patchPointI)
                    {
                        label pointI = subMeshPts[patchPointI];

                        typename Map<T>::const_iterator iter =
                            pointValues.find(pointI);

                        if (iter != pointValues.end())
                        {
                            subPatchInfo.insert(patchPointI, iter());
                        }
                    }

                    if (subPatchI != -1)
                    {
                        // Apply transform.
                        top
                        (
                            refCast<const coupledPolyPatch>
                            (
                                patches[subPatchI]
                            ),
                            subPatchInfo
                        );
                    }
                }

                OPstream toNeighb(Pstream::blocking, procPatch.neighbProcNo());
                toNeighb << patchInfo;
            }
        }


        // Receive and combine.

        forAll(patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].nPoints() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                List<Map<T> > nbrPatchInfo;
                {
                    IPstream fromNbr
                    (
                        Pstream::blocking,
                        procPatch.neighbProcNo()
                    );
                    fromNbr >> nbrPatchInfo;
                }

                forAll(procPatch.patchIDs(), subI)
                {
                    const labelList& subMeshPts =
                        procPatch.reverseSubMeshPoints()[subI];

                    const Map<T>& nbrSubInfo = nbrPatchInfo[subI];

                    forAllConstIter(typename Map<T>, nbrSubInfo, iter)
                    {
                        label meshPointI = subMeshPts[iter.key()];
                        combine
                        (
                            pointValues,
                            cop,
                            meshPointI,
                            iter()
                        );
                    }
                }
            }
        }
    }

    // Do the cyclics.
    forAll(patches, patchI)
    {
        if (isA<cyclicPolyPatch>(patches[patchI]))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchI]);

            if (cycPatch.owner())
            {
                // Owner does all.
                const cyclicPolyPatch& nbrPatch = cycPatch.neighbPatch();
                const edgeList& coupledPoints = cycPatch.coupledPoints();
                const labelList& meshPtsA = cycPatch.meshPoints();
                const labelList& meshPtsB = nbrPatch.meshPoints();

                // Extract local values. Create map from nbrPoint to value.
                Map<T> half0Values(meshPtsA.size() / 20);
                Map<T> half1Values(half0Values.size());

                forAll(coupledPoints, i)
                {
                    const edge& e = coupledPoints[i];

                    typename Map<T>::const_iterator point0Fnd =
                        pointValues.find(meshPtsA[e[0]]);

                    if (point0Fnd != pointValues.end())
                    {
                        half0Values.insert(i, point0Fnd());
                    }

                    typename Map<T>::const_iterator point1Fnd =
                        pointValues.find(meshPtsB[e[1]]);

                    if (point1Fnd != pointValues.end())
                    {
                        half1Values.insert(i, point1Fnd());
                    }
                }

                // Transform
                top(cycPatch, half0Values);
                top(nbrPatch, half1Values);

                forAll(coupledPoints, i)
                {
                    const edge& e = coupledPoints[i];

                    typename Map<T>::const_iterator half0Fnd =
                        half0Values.find(i);

                    if (half0Fnd != half0Values.end())
                    {
                        combine
                        (
                            pointValues,
                            cop,
                            meshPtsB[e[1]],
                            half0Fnd()
                        );
                    }

                    typename Map<T>::const_iterator half1Fnd =
                        half1Values.find(i);

                    if (half1Fnd != half1Values.end())
                    {
                        combine
                        (
                            pointValues,
                            cop,
                            meshPtsA[e[0]],
                            half1Fnd()
                        );
                    }
                }
            }
        }
    }

    // Synchronize multiple shared points.
    const globalMeshData& pd = mesh.globalData();

    if (pd.nGlobalPoints() > 0)
    {
        // meshPoint per local index
        const labelList& sharedPtLabels = pd.sharedPointLabels();
        // global shared index per local index
        const labelList& sharedPtAddr = pd.sharedPointAddr();

        // Values on shared points. Keyed on global shared index.
        Map<T> sharedPointValues(sharedPtAddr.size());


        // Fill my entries in the shared points
        forAll(sharedPtLabels, i)
        {
            label meshPointI = sharedPtLabels[i];

            typename Map<T>::const_iterator fnd =
                pointValues.find(meshPointI);

            if (fnd != pointValues.end())
            {
                combine
                (
                    sharedPointValues,
                    cop,
                    sharedPtAddr[i],    // index
                    fnd()               // value
                );
            }
        }

        // Reduce on master.

        if (Pstream::parRun())
        {
            if (Pstream::master())
            {
                // Receive the edges using shared points from the slave.
                for
                (
                    int slave=Pstream::firstSlave();
                    slave<=Pstream::lastSlave();
                    slave++
                )
                {
                    IPstream fromSlave(Pstream::blocking, slave);
                    Map<T> nbrValues(fromSlave);

                    // Merge neighbouring values with my values
                    forAllConstIter(typename Map<T>, nbrValues, iter)
                    {
                        combine
                        (
                            sharedPointValues,
                            cop,
                            iter.key(), // edge
                            iter()      // value
                        );
                    }
                }

                // Send back
                for
                (
                    int slave=Pstream::firstSlave();
                    slave<=Pstream::lastSlave();
                    slave++
                )
                {
                    OPstream toSlave(Pstream::blocking, slave);
                    toSlave << sharedPointValues;
                }
            }
            else
            {
                // Send to master
                {
                    OPstream toMaster(Pstream::blocking, Pstream::masterNo());
                    toMaster << sharedPointValues;
                }
                // Receive merged values
                {
                    IPstream fromMaster
                    (
                        Pstream::blocking,
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

        forAllConstIter(Map<label>, sharedToMeshPoint, iter)
        {
            // Do I have a value for my shared point
            typename Map<T>::const_iterator sharedFnd =
                sharedPointValues.find(iter.key());

            if (sharedFnd != sharedPointValues.end())
            {
                combine
                (
                    pointValues,
                    cop,
                    iter(),     // index
                    sharedFnd() // value
                );
            }
        }
    }
}


template <class T, class CombineOp, class TransformOp>
void Foam::syncTools::syncEdgeMap
(
    const polyMesh& mesh,
    EdgeMap<T>& edgeValues,
    const CombineOp& cop,
    const TransformOp& top
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (!hasCouples(patches))
    {
        return;
    }


    // Do synchronisation without constructing globalEdge addressing
    // (since this constructs mesh edge addressing)


    // Swap proc patch info
    // ~~~~~~~~~~~~~~~~~~~~

    if (Pstream::parRun())
    {
        // Send

        forAll(patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].nEdges() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                // Data per sub patch
                List<EdgeMap<T> > patchInfo(procPatch.patchIDs().size());

                forAll(procPatch.patchIDs(), subI)
                {
                    label subPatchI = procPatch.patchIDs()[subI];
                    label subStart = procPatch.starts()[subI];
                    label subSize = procPatch.starts()[subI+1]-subStart;

                    primitivePatch subPatch
                    (
                        faceSubList
                        (
                            procPatch,
                            subSize,
                            subStart
                        ),
                        procPatch.points()
                    );  
                    const labelList& subMeshPts = subPatch.meshPoints();
                    const edgeList& subEdges = subPatch.edges();

                    EdgeMap<T>& subPatchInfo = patchInfo[subI];
                    subPatchInfo.resize(subEdges.size());

                    forAll(subEdges, i)
                    {
                        const edge& e = subEdges[i];
                        const edge meshEdge(subMeshPts[e[0]], subMeshPts[e[1]]);

                        typename EdgeMap<T>::const_iterator iter =
                            edgeValues.find(meshEdge);

                        if (iter != edgeValues.end())
                        {
                            subPatchInfo.insert(e, iter());
                        }
                    }

                    if (subPatchI != -1)
                    {
                        // Apply transform.
                        top
                        (
                            refCast<const coupledPolyPatch>
                            (
                                patches[subPatchI]
                            ),
                            subPatchInfo
                        );
                    }
                }

                OPstream toNeighb(Pstream::blocking, procPatch.neighbProcNo());
                toNeighb << patchInfo;
            }
        }


        // Receive and combine.

        forAll(patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].nEdges() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                List<EdgeMap<T> > nbrPatchInfo;
                {
                    IPstream fromNbr
                    (
                        Pstream::blocking,
                        procPatch.neighbProcNo()
                    );
                    fromNbr >> nbrPatchInfo;
                }

                forAll(procPatch.patchIDs(), subI)
                {
                    const labelList& subMeshPts =
                        procPatch.subMeshPoints()[subI];

                    const EdgeMap<T>& nbrSubInfo = nbrPatchInfo[subI];

                    forAllConstIter(typename EdgeMap<T>, nbrSubInfo, iter)
                    {
                        const edge& e = iter.key();
                        const edge meshEdge(subMeshPts[e[0]], subMeshPts[e[1]]);

                        combine
                        (
                            edgeValues,
                            cop,
                            meshEdge,   // edge
                            iter()      // value
                        );
                    }
                }
            }
        }
    }


    // Swap cyclic info
    // ~~~~~~~~~~~~~~~~

    forAll(patches, patchI)
    {
        if (isA<cyclicPolyPatch>(patches[patchI]))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchI]);

            if (cycPatch.owner())
            {
                // Owner does all.

                const edgeList& coupledEdges = cycPatch.coupledEdges();
                const labelList& meshPtsA = cycPatch.meshPoints();
                const edgeList& edgesA = cycPatch.edges();
                const cyclicPolyPatch& nbrPatch = cycPatch.neighbPatch();
                const labelList& meshPtsB = nbrPatch.meshPoints();
                const edgeList& edgesB = nbrPatch.edges();

                // Extract local values. Create map from edge to value.
                Map<T> half0Values(edgesA.size() / 20);
                Map<T> half1Values(half0Values.size());

                forAll(coupledEdges, i)
                {
                    const edge& twoEdges = coupledEdges[i];

                    {
                        const edge& e0 = edgesA[twoEdges[0]];
                        const edge meshEdge0(meshPtsA[e0[0]], meshPtsA[e0[1]]);

                        typename EdgeMap<T>::const_iterator iter =
                            edgeValues.find(meshEdge0);

                        if (iter != edgeValues.end())
                        {
                            half0Values.insert(i, iter());
                        }
                    }
                    {
                        const edge& e1 = edgesB[twoEdges[1]];
                        const edge meshEdge1(meshPtsB[e1[0]], meshPtsB[e1[1]]);

                        typename EdgeMap<T>::const_iterator iter =
                            edgeValues.find(meshEdge1);

                        if (iter != edgeValues.end())
                        {
                            half1Values.insert(i, iter());
                        }
                    }
                }

                // Transform
                top(cycPatch, half0Values);
                top(nbrPatch, half1Values);


                // Extract and combine information

                forAll(coupledEdges, i)
                {
                    const edge& twoEdges = coupledEdges[i];

                    typename Map<T>::const_iterator half1Fnd =
                        half1Values.find(i);

                    if (half1Fnd != half1Values.end())
                    {
                        const edge& e0 = edgesA[twoEdges[0]];
                        const edge meshEdge0(meshPtsA[e0[0]], meshPtsA[e0[1]]);

                        combine
                        (
                            edgeValues,
                            cop,
                            meshEdge0,  // edge
                            half1Fnd()  // value
                        );
                    }

                    typename Map<T>::const_iterator half0Fnd =
                        half0Values.find(i);
                    if (half0Fnd != half0Values.end())
                    {
                        const edge& e1 = edgesB[twoEdges[1]];
                        const edge meshEdge1(meshPtsB[e1[0]], meshPtsB[e1[1]]);

                        combine
                        (
                            edgeValues,
                            cop,
                            meshEdge1,  // edge
                            half0Fnd()  // value
                        );
                    }
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
    for (label faceI = mesh.nInternalFaces(); faceI < mesh.nFaces(); faceI++)
    {
        const face& f = mesh.faces()[faceI];

        forAll(f, fp)
        {
            label v0 = f[fp];
            label v1 = f[f.fcIndex(fp)];

            Map<label>::const_iterator v0Fnd = meshToShared.find(v0);

            if (v0Fnd != meshToShared.end())
            {
                Map<label>::const_iterator v1Fnd = meshToShared.find(v1);

                if (v1Fnd != meshToShared.end())
                {
                    const edge meshEdge(v0, v1);

                    // edge in shared point labels
                    const edge sharedEdge(v0Fnd(), v1Fnd());

                    // Store mesh edge as a potential shared edge.
                    potentialSharedEdge.insert(sharedEdge, meshEdge);

                    typename EdgeMap<T>::const_iterator edgeFnd =
                        edgeValues.find(meshEdge);

                    if (edgeFnd != edgeValues.end())
                    {
                        // edge exists in edgeValues. See if already in map
                        // (since on same processor, e.g. cyclic)
                        combine
                        (
                            sharedEdgeValues,
                            cop,
                            sharedEdge, // edge
                            edgeFnd()   // value
                        );
                    }
                }
            }
        }
    }


    // Now sharedEdgeValues will contain per potential sharedEdge the value.
    // (potential since an edge having two shared points is not nessecary a
    //  shared edge).
    // Reduce this on the master.

    if (Pstream::parRun())
    {
        if (Pstream::master())
        {
            // Receive the edges using shared points from the slave.
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                IPstream fromSlave(Pstream::blocking, slave);
                EdgeMap<T> nbrValues(fromSlave);

                // Merge neighbouring values with my values
                forAllConstIter(typename EdgeMap<T>, nbrValues, iter)
                {
                    combine
                    (
                        sharedEdgeValues,
                        cop,
                        iter.key(), // edge
                        iter()      // value
                    );
                }
            }

            // Send back
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {

                OPstream toSlave(Pstream::blocking, slave);
                toSlave << sharedEdgeValues;
            }
        }
        else
        {
            // Send to master
            {
                OPstream toMaster(Pstream::blocking, Pstream::masterNo());
                toMaster << sharedEdgeValues;
            }
            // Receive merged values
            {
                IPstream fromMaster(Pstream::blocking, Pstream::masterNo());
                fromMaster >> sharedEdgeValues;
            }
        }
    }


    // Merge sharedEdgeValues (keyed on sharedPointAddr) into edgeValues
    // (keyed on mesh points).

    // Loop over all my shared edges.
    forAllConstIter(typename EdgeMap<edge>, potentialSharedEdge, iter)
    {
        const edge& sharedEdge = iter.key();
        const edge& meshEdge = iter();

        // Do I have a value for the shared edge?
        typename EdgeMap<T>::const_iterator sharedFnd =
            sharedEdgeValues.find(sharedEdge);

        if (sharedFnd != sharedEdgeValues.end())
        {
            combine
            (
                edgeValues,
                cop,
                meshEdge,       // edge
                sharedFnd()     // value
            );
        }
    }
}


template <class T, class CombineOp, class TransformOp>
void Foam::syncTools::syncPointList
(
    const polyMesh& mesh,
    UList<T>& pointValues,
    const CombineOp& cop,
    const T& nullValue,
    const TransformOp& top
)
{
    if (pointValues.size() != mesh.nPoints())
    {
        FatalErrorIn
        (
            "syncTools<class T, class CombineOp>::syncPointList"
            "(const polyMesh&, UList<T>&, const CombineOp&, const T&"
            ", const bool)"
        )   << "Number of values " << pointValues.size()
            << " is not equal to the number of points in the mesh "
            << mesh.nPoints() << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (!hasCouples(patches))
    {
        return;
    }

    if (Pstream::parRun())
    {
        // Send

        forAll(patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].nPoints() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                //TBD: optimisation for single subpatch.

                // Data per sub patch
                List<List<T> > patchInfo(procPatch.patchIDs().size());

                forAll(procPatch.patchIDs(), subI)
                {
                    label subPatchI = procPatch.patchIDs()[subI];

                    const labelList& subMeshPts =
                        procPatch.subMeshPoints()[subI];

                    List<T>& subPatchInfo = patchInfo[subI];
                    subPatchInfo.setSize(subMeshPts.size());
                    forAll(subPatchInfo, i)
                    {
                        subPatchInfo[i] = pointValues[subMeshPts[i]];
                    }

                    if (subPatchI != -1)
                    {
                        // Apply transform.
                        SubField<T> slice(subPatchInfo, subPatchInfo.size());
                        top
                        (
                            refCast<const coupledPolyPatch>
                            (
                                patches[subPatchI]
                            ),
                            reinterpret_cast<Field<T>&>(slice)
                        );
                    }
                }

                OPstream toNbr(Pstream::blocking, procPatch.neighbProcNo());
                toNbr << patchInfo;
            }
        }


        // Receive and combine.

        forAll(patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].nPoints() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                List<List<T> > nbrPatchInfo;
                {
                    // We do not know the number of points on the other side
                    // so cannot use Pstream::read.
                    IPstream fromNbr
                    (
                        Pstream::blocking,
                        procPatch.neighbProcNo()
                    );
                    fromNbr >> nbrPatchInfo;
                }

                forAll(procPatch.patchIDs(), subI)
                {
                    const labelList& subMeshPts =
                        procPatch.reverseSubMeshPoints()[subI];
                    const List<T>& nbrSubInfo = nbrPatchInfo[subI];

                    forAll(subMeshPts, i)
                    {
                        label meshPointI = subMeshPts[i];
                        cop(pointValues[meshPointI], nbrSubInfo[i]);
                    }
                }
            }
        }
    }

    // Do the cyclics.
    forAll(patches, patchI)
    {
        if (isA<cyclicPolyPatch>(patches[patchI]))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchI]);

            if (cycPatch.owner())
            {
                // Owner does all.

                const edgeList& coupledPoints = cycPatch.coupledPoints();
                const labelList& meshPts = cycPatch.meshPoints();
                const cyclicPolyPatch& nbrPatch = cycPatch.neighbPatch();
                const labelList& nbrMeshPoints = nbrPatch.meshPoints();

                Field<T> half0Values(coupledPoints.size());
                Field<T> half1Values(coupledPoints.size());

                forAll(coupledPoints, i)
                {
                    const edge& e = coupledPoints[i];
                    half0Values[i] = pointValues[meshPts[e[0]]];
                    half1Values[i] = pointValues[nbrMeshPoints[e[1]]];
                }

                //SubField<T> slice0(half0Values, half0Values.size());
                //SubField<T> slice1(half1Values, half1Values.size());
                //top(cycPatch, reinterpret_cast<Field<T>&>(slice0));
                //top(nbrPatch, reinterpret_cast<Field<T>&>(slice1));

                top(cycPatch, half0Values);
                top(nbrPatch, half1Values);

                forAll(coupledPoints, i)
                {
                    const edge& e = coupledPoints[i];
                    cop(pointValues[meshPts[e[0]]], half1Values[i]);
                    cop(pointValues[nbrMeshPoints[e[1]]], half0Values[i]);
                }
            }
        }
    }

    // Synchronize multiple shared points.
    const globalMeshData& pd = mesh.globalData();

    if (pd.nGlobalPoints() > 0)
    {
        // Values on shared points.
        List<T> sharedPts(pd.nGlobalPoints(), nullValue);

        forAll(pd.sharedPointLabels(), i)
        {
            label meshPointI = pd.sharedPointLabels()[i];
            // Fill my entries in the shared points
            sharedPts[pd.sharedPointAddr()[i]] = pointValues[meshPointI];
        }

        // Combine on master.
        Pstream::listCombineGather(sharedPts, cop);
        Pstream::listCombineScatter(sharedPts);

        // Now we will all have the same information. Merge it back with
        // my local information.
        forAll(pd.sharedPointLabels(), i)
        {
            label meshPointI = pd.sharedPointLabels()[i];
            pointValues[meshPointI] = sharedPts[pd.sharedPointAddr()[i]];
        }
    }
}


template <class T, class CombineOp, class TransformOp>
void Foam::syncTools::syncPointList
(
    const polyMesh& mesh,
    const labelList& meshPoints,
    UList<T>& pointValues,
    const CombineOp& cop,
    const T& nullValue,
    const TransformOp& top
)
{
    if (pointValues.size() != meshPoints.size())
    {
        FatalErrorIn
        (
            "syncTools<class T, class CombineOp>::syncPointList"
            "(const polyMesh&, const labelList&, UList<T>&, const CombineOp&"
            ", const T&, const bool)"
        )   << "Number of values " << pointValues.size()
            << " is not equal to the number of points "
            << meshPoints.size() << abort(FatalError);
    }

    if (!hasCouples(mesh.boundaryMesh()))
    {
        return;
    }

    List<T> meshValues(mesh.nPoints(), nullValue);

    forAll(meshPoints, i)
    {
        meshValues[meshPoints[i]] = pointValues[i];
    }

    syncTools::syncPointList
    (
        mesh,
        meshValues,
        cop,            // combine op
        nullValue,      // null value
        top             // position or field
    );

    forAll(meshPoints, i)
    {
        pointValues[i] = meshValues[meshPoints[i]];
    }
}


template <class T, class CombineOp, class TransformOp>
void Foam::syncTools::syncEdgeList
(
    const polyMesh& mesh,
    UList<T>& edgeValues,
    const CombineOp& cop,
    const T& nullValue,
    const TransformOp& top
)
{
    if (edgeValues.size() != mesh.nEdges())
    {
        FatalErrorIn
        (
            "syncTools<class T, class CombineOp>::syncEdgeList"
            "(const polyMesh&, UList<T>&, const CombineOp&, const T&"
            ", const bool)"
        )   << "Number of values " << edgeValues.size()
            << " is not equal to the number of edges in the mesh "
            << mesh.nEdges() << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (!hasCouples(patches))
    {
        return;
    }

    if (Pstream::parRun())
    {
        // Send

        forAll(patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].nEdges() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                //TBD: optimisation for single subpatch.

                // Data per sub patch
                List<List<T> > patchInfo(procPatch.patchIDs().size());

                forAll(procPatch.patchIDs(), subI)
                {
                    label subPatchI = procPatch.patchIDs()[subI];
                    label subStart = procPatch.starts()[subI];
                    label subSize = procPatch.starts()[subI+1]-subStart;

                    primitivePatch subPatch
                    (
                        faceSubList
                        (
                            procPatch,
                            subSize,
                            subStart
                        ),
                        procPatch.points()
                    );
                    labelList subMeshEdges
                    (
                        subPatch.meshEdges
                        (
                            mesh.edges(),
                            mesh.pointEdges()
                        )
                    );

                    List<T>& subPatchInfo = patchInfo[subI];
                    subPatchInfo.setSize(subMeshEdges.size());
                    forAll(subPatchInfo, i)
                    {
                        subPatchInfo[i] = edgeValues[subMeshEdges[i]];
                    }

                    if (subPatchI != -1)
                    {
                        // Apply transform.
                        SubField<T> slice(subPatchInfo, subPatchInfo.size());
                        top
                        (
                            refCast<const coupledPolyPatch>
                            (
                                patches[subPatchI]
                            ),
                            reinterpret_cast<Field<T>&>(slice)
                        );
                    }
                }

                OPstream toNbr(Pstream::blocking, procPatch.neighbProcNo());
                toNbr << patchInfo;
           }
        }

        // Receive and combine.

        forAll(patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].nEdges() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                // Receive from neighbour. 
                List<List<T> > nbrPatchInfo(procPatch.nEdges());

                {
                    IPstream fromNeighb
                    (
                        Pstream::blocking,
                        procPatch.neighbProcNo()
                    );
                    fromNeighb >> nbrPatchInfo;
                }

                forAll(procPatch.patchIDs(), subI)
                {
                    label subStart = procPatch.starts()[subI];
                    label subSize = procPatch.starts()[subI+1]-subStart;

                    faceList reverseFaces(subSize);
                    forAll(reverseFaces, i)
                    {
                        reverseFaces[i] = procPatch[subStart+i].reverseFace();
                    }
                    primitivePatch subPatch
                    (
                        faceSubList
                        (
                            reverseFaces,
                            reverseFaces.size()
                        ),
                        procPatch.points()
                    );
                    labelList subMeshEdges
                    (
                        subPatch.meshEdges
                        (
                            mesh.edges(),
                            mesh.pointEdges()
                        )
                    );

                    const List<T>& nbrSubInfo = nbrPatchInfo[subI];

                    forAll(subMeshEdges, i)
                    {
                        label meshEdgeI = subMeshEdges[i];
                        cop(edgeValues[meshEdgeI], nbrSubInfo[i]);
                    }
                }
            }
        }
    }

    // Do the cyclics.
    forAll(patches, patchI)
    {
        if (isA<cyclicPolyPatch>(patches[patchI]))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchI]);

            if (cycPatch.owner())
            {
                // Owner does all.
                const edgeList& coupledEdges = cycPatch.coupledEdges();
                const labelList& meshEdges = cycPatch.meshEdges();
                const cyclicPolyPatch& nbrPatch = cycPatch.neighbPatch();
                const labelList& nbrMeshEdges = nbrPatch.meshEdges();

                Field<T> half0Values(coupledEdges.size());
                Field<T> half1Values(coupledEdges.size());

                forAll(coupledEdges, i)
                {
                    const edge& e = coupledEdges[i];
                    half0Values[i] = edgeValues[meshEdges[e[0]]];
                    half1Values[i] = edgeValues[nbrMeshEdges[e[1]]];
                }

                //SubField<T> slice0(half0Values, half0Values.size());
                //SubField<T> slice1(half1Values, half1Values.size());
                //top(cycPatch, reinterpret_cast<Field<T>&>(slice0));
                //top(nbrPatch, reinterpret_cast<Field<T>&>(slice1));

                top(cycPatch, half0Values);
                top(nbrPatch, half1Values);

                forAll(coupledEdges, i)
                {
                    const edge& e = coupledEdges[i];
                    cop(edgeValues[meshEdges[e[0]]], half1Values[i]);
                    cop(edgeValues[nbrMeshEdges[e[1]]], half0Values[i]);
                }
            }
        }
    }

    //- Note: hasTransformation is only used for warning messages so
    //  reduction not strictly nessecary.
    //reduce(hasTransformation, orOp<bool>());

    // Do the multiple shared edges
    const globalMeshData& pd = mesh.globalData();

    if (pd.nGlobalEdges() > 0)
    {
        // Values on shared edges.
        List<T> sharedPts(pd.nGlobalEdges(), nullValue);

        forAll(pd.sharedEdgeLabels(), i)
        {
            label meshEdgeI = pd.sharedEdgeLabels()[i];

            // Fill my entries in the shared edges
            sharedPts[pd.sharedEdgeAddr()[i]] = edgeValues[meshEdgeI];
        }

        // Combine on master.
        Pstream::listCombineGather(sharedPts, cop);
        Pstream::listCombineScatter(sharedPts);

        // Now we will all have the same information. Merge it back with
        // my local information.
        forAll(pd.sharedEdgeLabels(), i)
        {
            label meshEdgeI = pd.sharedEdgeLabels()[i];
            edgeValues[meshEdgeI] = sharedPts[pd.sharedEdgeAddr()[i]];
        }
    }
}


template <class T, class CombineOp, class TransformOp>
void Foam::syncTools::syncBoundaryFaceList
(
    const polyMesh& mesh,
    UList<T>& faceValues,
    const CombineOp& cop,
    const TransformOp& top
)
{
    const label nBFaces = mesh.nFaces() - mesh.nInternalFaces();

    if (faceValues.size() != nBFaces)
    {
        FatalErrorIn
        (
            "syncTools<class T, class CombineOp>::syncBoundaryFaceList"
            "(const polyMesh&, UList<T>&, const CombineOp&"
            ", const bool)"
        )   << "Number of values " << faceValues.size()
            << " is not equal to the number of boundary faces in the mesh "
            << nBFaces << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (!hasCouples(patches))
    {
        return;
    }


    if (Pstream::parRun())
    {
        // Send

        forAll(patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].size() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                label patchStart = procPatch.start()-mesh.nInternalFaces();

                if (contiguous<T>())
                {
                    OPstream::write
                    (
                        Pstream::blocking,
                        procPatch.neighbProcNo(),
                        reinterpret_cast<const char*>(&faceValues[patchStart]),
                        procPatch.size()*sizeof(T)
                    );
                }
                else
                {
                    OPstream toNbr(Pstream::blocking, procPatch.neighbProcNo());
                    toNbr <<
                        SubList<T>(faceValues, procPatch.size(), patchStart);
                }
            }
        }


        // Receive and combine.

        forAll(patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].size() > 0
            )
            {
                const processorPolyPatch& ppp =
                    refCast<const processorPolyPatch>(patches[patchI]);

                List<T> nbrPatchInfo(ppp.size());

                if (contiguous<T>())
                {
                    IPstream::read
                    (
                        Pstream::blocking,
                        ppp.neighbProcNo(),
                        reinterpret_cast<char*>(nbrPatchInfo.begin()),
                        nbrPatchInfo.byteSize()
                    );
                }
                else
                {
                    IPstream fromNeighb
                    (
                        Pstream::blocking,
                        ppp.neighbProcNo()
                    );
                    fromNeighb >> nbrPatchInfo;
                }

                // Apply any transforms
                forAll(ppp.patchIDs(), i)
                {
                    label subPatchI = ppp.patchIDs()[i];
                    label subStart = ppp.starts()[i];
                    label subSize = ppp.starts()[i+1]-subStart;

                    if (subPatchI != -1)
                    {
                        SubField<T> slice(nbrPatchInfo, subStart, subSize);
                        top
                        (
                            refCast<const coupledPolyPatch>
                            (
                                patches[subPatchI]
                            ),
                            reinterpret_cast<Field<T>&>(slice)
                        );
                    }
                }

                label bFaceI = ppp.start()-mesh.nInternalFaces();

                forAll(nbrPatchInfo, i)
                {
                    cop(faceValues[bFaceI++], nbrPatchInfo[i]);
                }
            }
        }
    }

    // Do the cyclics.
    forAll(patches, patchI)
    {
        if (isA<cyclicPolyPatch>(patches[patchI]))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchI]);

            if (cycPatch.owner())
            {
                // Owner does all.
                const cyclicPolyPatch& nbrPatch = cycPatch.neighbPatch();
                label ownStart = cycPatch.start()-mesh.nInternalFaces();
                label nbrStart = nbrPatch.start()-mesh.nInternalFaces();

                label sz = cycPatch.size();

                // Transform (copy of) data on both sides
                Field<T> ownVals(SubList<T>(faceValues, sz, ownStart));
                //SubField<T> ownSlice(ownVals, ownVals.size());
                top(cycPatch, ownVals);

                Field<T> nbrVals(SubList<T>(faceValues, sz, nbrStart));
                //SubField<T> nbrSlice(nbrVals, nbrVals.size());
                top(nbrPatch, nbrVals);

                label i0 = ownStart;
                forAll(nbrVals, i)
                {
                    cop(faceValues[i0++], nbrVals[i]);
                }

                label i1 = nbrStart;
                forAll(ownVals, i)
                {
                    cop(faceValues[i1++], ownVals[i]);
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <unsigned nBits, class CombineOp>
void Foam::syncTools::syncFaceList
(
    const polyMesh& mesh,
    PackedList<nBits>& faceValues,
    const CombineOp& cop
)
{
    if (faceValues.size() != mesh.nFaces())
    {
        FatalErrorIn
        (
            "syncTools<unsigned nBits, class CombineOp>::syncFaceList"
            "(const polyMesh&, PackedList<nBits>&, const CombineOp&)"
        )   << "Number of values " << faceValues.size()
            << " is not equal to the number of faces in the mesh "
            << mesh.nFaces() << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (!hasCouples(patches))
    {
        return;
    }

    // Patch data (proc patches only) for ease of sending/receiving.
    // Note:should just send slice of packedList itself.
    List<List<unsigned int> > patchValues(patches.size());

    if (Pstream::parRun())
    {
        // Send

        forAll(patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].size() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                patchValues[patchI].setSize(procPatch.size());
                forAll(procPatch, i)
                {
                    patchValues[patchI][i] = faceValues[procPatch.start()+i];
                }

                OPstream toNbr(Pstream::blocking, procPatch.neighbProcNo());
                toNbr << patchValues[patchI];
            }
        }


        // Receive and combine.

        forAll(patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].size() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                {
                    IPstream fromNbr
                    (
                        Pstream::blocking,
                        procPatch.neighbProcNo()
                    );
                    fromNbr >> patchValues[patchI];
                }

                // Combine (bitwise)
                forAll(procPatch, i)
                {
                    unsigned int patchVal = patchValues[patchI][i];
                    label meshFaceI = procPatch.start()+i;
                    unsigned int faceVal = faceValues[meshFaceI];
                    cop(faceVal, patchVal);
                    faceValues[meshFaceI] = faceVal;
                }
            }
        }
    }

    // Do the cyclics.
    forAll(patches, patchI)
    {
        if (isA<cyclicPolyPatch>(patches[patchI]))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchI]);

            if (cycPatch.owner())
            {
                // Owner does all.
                const cyclicPolyPatch& nbrPatch = cycPatch.neighbPatch();

                for (label i = 0; i < cycPatch.size(); i++)
                {
                    label meshFace0 = cycPatch.start()+i;
                    unsigned int val0 = faceValues[meshFace0];
                    label meshFace1 = nbrPatch.start()+i;
                    unsigned int val1 = faceValues[meshFace1];

                    unsigned int t = val0;
                    cop(t, val1);
                    faceValues[meshFace0] = t;

                    cop(val1, val0);                
                    faceValues[meshFace1] = val1;
                }
            }
        }
    }
}


template <unsigned nBits>
void Foam::syncTools::swapFaceList
(
    const polyMesh& mesh,
    PackedList<nBits>& faceValues
)
{
    syncFaceList(mesh, faceValues, eqOp<unsigned int>());
}


template <unsigned nBits, class CombineOp>
void Foam::syncTools::syncPointList
(
    const polyMesh& mesh,
    PackedList<nBits>& pointValues,
    const CombineOp& cop,
    const unsigned int nullValue
)
{
    if (pointValues.size() != mesh.nPoints())
    {
        FatalErrorIn
        (
            "syncTools<unsigned nBits, class CombineOp>::syncPointList"
            "(const polyMesh&, PackedList<nBits>&, const CombineOp&"
            ", const unsigned int)"
        )   << "Number of values " << pointValues.size()
            << " is not equal to the number of points in the mesh "
            << mesh.nPoints() << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (!hasCouples(patches))
    {
        return;
    }

    if (Pstream::parRun())
    {
        // Send

        forAll(patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].nPoints() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                List<List<unsigned int> > patchInfo(procPatch.patchIDs().size());

                forAll(procPatch.patchIDs(), subI)
                {
                    const labelList& subMeshPts =
                        procPatch.subMeshPoints()[subI];

                    List<unsigned int>& subInfo = patchInfo[subI];
                    subInfo.setSize(subMeshPts.size());
                    forAll(subInfo, i)
                    {
                        subInfo[i] = pointValues[subMeshPts[i]];
                    }
                }

                OPstream toNbr(Pstream::blocking, procPatch.neighbProcNo());
                toNbr << patchInfo;
            }
        }


        // Receive and combine.

        forAll(patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].nPoints() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                List<List<unsigned int> > nbrPatchInfo(procPatch.nPoints());
                {
                    // We do not know the number of points on the other side
                    // so cannot use Pstream::read.
                    IPstream fromNbr
                    (
                        Pstream::blocking,
                        procPatch.neighbProcNo()
                    );
                    fromNbr >> nbrPatchInfo;
                }

                forAll(procPatch.patchIDs(), subI)
                {
                    const labelList& subMeshPts =
                        procPatch.reverseSubMeshPoints()[subI];
                    const List<unsigned int>& nbrSubInfo = nbrPatchInfo[subI];

                    forAll(subMeshPts, i)
                    {
                        label meshPointI = subMeshPts[i];
                        unsigned int pointVal = pointValues[meshPointI];
                        cop(pointVal, nbrSubInfo[i]);
                        pointValues[meshPointI] = pointVal;
                    }
                }
            }
        }
    }

    // Do the cyclics.
    forAll(patches, patchI)
    {
        if (isA<cyclicPolyPatch>(patches[patchI]))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchI]);

            if (cycPatch.owner())
            {
                // Owner does all.

                const edgeList& coupledPoints = cycPatch.coupledPoints();
                const labelList& meshPts = cycPatch.meshPoints();
                const cyclicPolyPatch& nbrPatch = cycPatch.neighbPatch();
                const labelList& nbrMeshPts = nbrPatch.meshPoints();

                forAll(coupledPoints, i)
                {
                    const edge& e = coupledPoints[i];
                    label meshPoint0 = meshPts[e[0]];
                    label meshPoint1 = nbrMeshPts[e[1]];

                    unsigned int val0 = pointValues[meshPoint0];
                    unsigned int val1 = pointValues[meshPoint1];
                    unsigned int t = val0;

                    cop(val0, val1);
                    pointValues[meshPoint0] = val0;
                    cop(val1, t);
                    pointValues[meshPoint1] = val1;
                }
            }
        }
    }

    // Synchronize multiple shared points.
    const globalMeshData& pd = mesh.globalData();

    if (pd.nGlobalPoints() > 0)
    {
        // Values on shared points. Use unpacked storage for ease!
        List<unsigned int> sharedPts(pd.nGlobalPoints(), nullValue);

        forAll(pd.sharedPointLabels(), i)
        {
            label meshPointI = pd.sharedPointLabels()[i];
            // Fill my entries in the shared points
            sharedPts[pd.sharedPointAddr()[i]] = pointValues[meshPointI];
        }

        // Combine on master.
        Pstream::listCombineGather(sharedPts, cop);
        Pstream::listCombineScatter(sharedPts);

        // Now we will all have the same information. Merge it back with
        // my local information.
        forAll(pd.sharedPointLabels(), i)
        {
            label meshPointI = pd.sharedPointLabels()[i];
            pointValues[meshPointI] = sharedPts[pd.sharedPointAddr()[i]];
        }
    }
}


template <unsigned nBits, class CombineOp>
void Foam::syncTools::syncEdgeList
(
    const polyMesh& mesh,
    PackedList<nBits>& edgeValues,
    const CombineOp& cop,
    const unsigned int nullValue
)
{
    if (edgeValues.size() != mesh.nEdges())
    {
        FatalErrorIn
        (
            "syncTools<unsigned nBits, class CombineOp>::syncEdgeList"
            "(const polyMesh&, PackedList<nBits>&, const CombineOp&"
            ", const unsigned int)"
        )   << "Number of values " << edgeValues.size()
            << " is not equal to the number of edges in the mesh "
            << mesh.nEdges() << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (!hasCouples(patches))
    {
        return;
    }

    if (Pstream::parRun())
    {
        // Send

        forAll(patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].nEdges() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                List<List<unsigned int> > patchInfo(procPatch.patchIDs().size());

                forAll(procPatch.patchIDs(), subI)
                {
                    label subStart = procPatch.starts()[subI];
                    label subSize = procPatch.starts()[subI+1]-subStart;

                    primitivePatch subPatch
                    (
                        faceSubList
                        (
                            procPatch,
                            subSize,
                            subStart
                        ),
                        procPatch.points()
                    );
                    labelList subMeshEdges
                    (
                        subPatch.meshEdges
                        (
                            mesh.edges(),
                            mesh.pointEdges()
                        )
                    );

                    List<unsigned int>& subPatchInfo = patchInfo[subI];
                    subPatchInfo.setSize(subMeshEdges.size());
                    forAll(subPatchInfo, i)
                    {
                        subPatchInfo[i] = edgeValues[subMeshEdges[i]];
                    }
                }

                OPstream toNbr(Pstream::blocking, procPatch.neighbProcNo());
                toNbr << patchInfo;
            }
        }


        // Receive and combine.

        forAll(patches, patchI)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchI])
             && patches[patchI].nEdges() > 0
            )
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                // Receive from neighbour. 
                List<List<unsigned int> > nbrPatchInfo(procPatch.nEdges());

                {
                    IPstream fromNeighb
                    (
                        Pstream::blocking,
                        procPatch.neighbProcNo()
                    );
                    fromNeighb >> nbrPatchInfo;
                }

                forAll(procPatch.patchIDs(), subI)
                {
                    label subStart = procPatch.starts()[subI];
                    label subSize = procPatch.starts()[subI+1]-subStart;

                    faceList reverseFaces(subSize);
                    forAll(reverseFaces, i)
                    {
                        reverseFaces[i] = procPatch[subStart+i].reverseFace();
                    }
                    primitivePatch subPatch
                    (
                        faceSubList
                        (
                            reverseFaces,
                            reverseFaces.size()
                        ),
                        procPatch.points()
                    );
                    labelList subMeshEdges
                    (
                        subPatch.meshEdges
                        (
                            mesh.edges(),
                            mesh.pointEdges()
                        )
                    );

                    const List<unsigned int>& nbrSubInfo = nbrPatchInfo[subI];

                    forAll(subMeshEdges, i)
                    {
                        label meshEdgeI = subMeshEdges[i];
                        unsigned int edgeVal = edgeValues[meshEdgeI];
                        cop(edgeVal, nbrSubInfo[i]);
                        edgeValues[meshEdgeI] = edgeVal;
                    }
                }
            }
        }
    }

    // Do the cyclics.
    forAll(patches, patchI)
    {
        if (isA<cyclicPolyPatch>(patches[patchI]))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patches[patchI]);

            if (cycPatch.owner())
            {
                // Owner does all.
                const edgeList& coupledEdges = cycPatch.coupledEdges();
                const labelList& meshEdges = cycPatch.meshEdges();
                const cyclicPolyPatch& nbrPatch = cycPatch.neighbPatch();
                const labelList& nbrMeshEdges = nbrPatch.meshEdges();

                forAll(coupledEdges, i)
                {
                    const edge& e = coupledEdges[i];

                    label edge0 = meshEdges[e[0]];
                    label edge1 = nbrMeshEdges[e[1]];

                    unsigned int val0 = edgeValues[edge0];
                    unsigned int t = val0;
                    unsigned int val1 = edgeValues[edge1];

                    cop(t, val1);
                    edgeValues[edge0] = t;
                    cop(val1, val0);
                    edgeValues[edge1] = val1;
                }
            }
        }
    }

    // Synchronize multiple shared edges.
    const globalMeshData& pd = mesh.globalData();

    if (pd.nGlobalEdges() > 0)
    {
        // Values on shared edges. Use unpacked storage for ease!
        List<unsigned int> sharedPts(pd.nGlobalEdges(), nullValue);

        forAll(pd.sharedEdgeLabels(), i)
        {
            label meshEdgeI = pd.sharedEdgeLabels()[i];
            // Fill my entries in the shared edges
            sharedPts[pd.sharedEdgeAddr()[i]] = edgeValues[meshEdgeI];
        }

        // Combine on master.
        Pstream::listCombineGather(sharedPts, cop);
        Pstream::listCombineScatter(sharedPts);

        // Now we will all have the same information. Merge it back with
        // my local information.
        forAll(pd.sharedEdgeLabels(), i)
        {
            label meshEdgeI = pd.sharedEdgeLabels()[i];
            edgeValues[meshEdgeI] = sharedPts[pd.sharedEdgeAddr()[i]];
        }
    }
}


// ************************************************************************* //
