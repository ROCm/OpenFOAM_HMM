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

\*---------------------------------------------------------------------------*/

#include "globalPoints.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::globalPoints, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Total number of points on processor patches. Is upper limit for number
// of shared points
Foam::label Foam::globalPoints::countPatchPoints
(
    const polyBoundaryMesh& patches
)
{
    label nTotPoints = 0;

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if
        (
            (Pstream::parRun() && isA<processorPolyPatch>(pp))
         || isA<cyclicPolyPatch>(pp)
        )
        {
            nTotPoints += pp.nPoints();
        }
    }

    return nTotPoints;
}


// Collect all topological information about a point on a patch.
// (this information is the patch faces using the point and the relative
// position of the point in the face)
void Foam::globalPoints::addToSend
(
    const primitivePatch& pp,
    const label patchPointI,
    const labelList& knownInfo,

    DynamicList<label>& patchFaces,
    DynamicList<label>& indexInFace,
    DynamicList<labelList>& allInfo
)
{
    label meshPointI = pp.meshPoints()[patchPointI];

    // Add all faces using the point so we are sure we find it on the
    // other side.
    const labelList& pFaces = pp.pointFaces()[patchPointI];

    forAll(pFaces, i)
    {
        label patchFaceI = pFaces[i];

        const face& f = pp[patchFaceI];

        patchFaces.append(patchFaceI);
        indexInFace.append(findIndex(f, meshPointI));
        allInfo.append(knownInfo);
    }
}


// Add nbrInfo to myInfo. Return true if anything changed.
// nbrInfo is for a point a list of all the global points using it
bool Foam::globalPoints::mergeInfo
(
    const labelList& nbrInfo,
    labelList& myInfo
)
{
    labelList newInfo(myInfo);
    label newI = newInfo.size();
    newInfo.setSize(newI + nbrInfo.size());

    forAll(nbrInfo, i)
    {
        label index = findIndex(myInfo, nbrInfo[i]);

        if (index == -1)
        {
            newInfo[newI++] = nbrInfo[i];
        }
    }

    newInfo.setSize(newI);

    // Did anything change?
    bool anyChanged = (newI > myInfo.size());

    myInfo.transfer(newInfo);

    return anyChanged;
}


Foam::label Foam::globalPoints::meshToLocalPoint
(
    const Map<label>& meshToPatchPoint, // from mesh point to local numbering
    const label meshPointI
)
{
    return
    (
        meshToPatchPoint.size() == 0
      ? meshPointI
      : meshToPatchPoint[meshPointI]
    );
}


Foam::label Foam::globalPoints::localToMeshPoint
(
    const labelList& patchToMeshPoint,
    const label localPointI
)
{
    return
    (
        patchToMeshPoint.size() == 0
      ? localPointI
      : patchToMeshPoint[localPointI]
    );
}


// Updates database of current information on meshpoints with nbrInfo.
// Uses mergeInfo above. Returns true if data kept for meshPointI changed.
bool Foam::globalPoints::storeInfo
(
    const labelList& nbrInfo,
    const label localPointI
)
{
    label infoChanged = false;

    // Get the index into the procPoints list.
    Map<label>::iterator iter = meshToProcPoint_.find(localPointI);

    if (iter != meshToProcPoint_.end())
    {
        if (mergeInfo(nbrInfo, procPoints_[iter()]))
        {
            infoChanged = true;
        }
    }
    else
    {
        labelList knownInfo(1, globalIndices_.toGlobal(localPointI));

        if (mergeInfo(nbrInfo, knownInfo))
        {
            // Update addressing from into procPoints
            meshToProcPoint_.insert(localPointI, procPoints_.size());
            // Insert into list of equivalences.
            procPoints_.append(knownInfo);

            infoChanged = true;
        }
    }
    return infoChanged;
}


// Insert my own points into structure and mark as changed.
void Foam::globalPoints::initOwnPoints
(
    const Map<label>& meshToPatchPoint,
    const bool allPoints,
    labelHashSet& changedPoints
)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if
        (
            (Pstream::parRun() && isA<processorPolyPatch>(pp))
         || isA<cyclicPolyPatch>(pp)
        )
        {
            const labelList& meshPoints = pp.meshPoints();

            if (allPoints)
            {
                // All points on patch
                forAll(meshPoints, i)
                {
                    label meshPointI = meshPoints[i];
                    label localPointI = meshToLocalPoint
                    (
                        meshToPatchPoint,
                        meshPointI
                    );
                    labelList knownInfo
                    (
                        1,
                        globalIndices_.toGlobal(localPointI)
                    );

                    // Update addressing from point to index in procPoints
                    meshToProcPoint_.insert(localPointI, procPoints_.size());
                    // Insert into list of equivalences.
                    procPoints_.append(knownInfo);

                    // Update changedpoints info.
                    changedPoints.insert(localPointI);
                }
            }
            else
            {
                // Boundary points only
                const labelList& boundaryPoints = pp.boundaryPoints();

                forAll(boundaryPoints, i)
                {
                    label meshPointI = meshPoints[boundaryPoints[i]];
                    label localPointI = meshToLocalPoint
                    (
                        meshToPatchPoint,
                        meshPointI
                    );

                    labelList knownInfo
                    (
                        1,
                        globalIndices_.toGlobal(localPointI)
                    );

                    // Update addressing from point to index in procPoints
                    meshToProcPoint_.insert(localPointI, procPoints_.size());
                    // Insert into list of equivalences.
                    procPoints_.append(knownInfo);

                    // Update changedpoints info.
                    changedPoints.insert(localPointI);
                }
            }
        }
    }
}


// Send all my info on changedPoints_ to my neighbours.
void Foam::globalPoints::sendPatchPoints
(
    const Map<label>& meshToPatchPoint,
    PstreamBuffers& pBufs,
    const labelHashSet& changedPoints
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (Pstream::parRun() && isA<processorPolyPatch>(pp))
        {
            // Information to send:
            // patch face
            DynamicList<label> patchFaces(pp.nPoints());
            // index in patch face
            DynamicList<label> indexInFace(pp.nPoints());
            // all information I currently hold about this patchPoint
            DynamicList<labelList> allInfo(pp.nPoints());


            // Now collect information on all points mentioned in
            // changedPoints. Note that these points only should occur on
            // processorPatches (or rather this is a limitation!).

            const labelList& meshPoints = pp.meshPoints();

            forAll(meshPoints, patchPointI)
            {
                label meshPointI = meshPoints[patchPointI];
                label localPointI = meshToLocalPoint
                (
                    meshToPatchPoint,
                    meshPointI
                );

                if (changedPoints.found(localPointI))
                {
                    label index = meshToProcPoint_[localPointI];

                    const labelList& knownInfo = procPoints_[index];

                    // Add my information about localPointI to the send buffers
                    addToSend
                    (
                        pp,
                        patchPointI,
                        knownInfo,

                        patchFaces,
                        indexInFace,
                        allInfo
                    );
                }
            }

            // Send to neighbour
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(pp);

                if (debug)
                {
                    Pout<< " Sending to "
                        << procPatch.neighbProcNo() << "   point information:"
                        << patchFaces.size() << endl;
                }

                UOPstream toNeighbour(procPatch.neighbProcNo(), pBufs);
                toNeighbour << patchFaces << indexInFace << allInfo;
            }
        }
    }
}


// Receive all my neighbours' information and merge with mine.
// After finishing will have updated
// - procPoints_ : all neighbour information merged in.
// - meshToProcPoint_
// - changedPoints: all points for which something changed.
void Foam::globalPoints::receivePatchPoints
(
    const Map<label>& meshToPatchPoint,
    PstreamBuffers& pBufs,
    labelHashSet& changedPoints
)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Reset changed points
    changedPoints.clear();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (Pstream::parRun() && isA<processorPolyPatch>(pp))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pp);

            labelList patchFaces;
            labelList indexInFace;
            List<labelList> nbrInfo;

            {
                UIPstream fromNeighbour(procPatch.neighbProcNo(), pBufs);
                fromNeighbour >> patchFaces >> indexInFace >> nbrInfo;
            }

            if (debug)
            {
                Pout<< " Received from "
                    << procPatch.neighbProcNo() << "   point information:"
                    << patchFaces.size() << endl;
            }

            forAll(patchFaces, i)
            {
                const face& f = pp[patchFaces[i]];

                // Get index in this face from index on face on other side.
                label index = (f.size() - indexInFace[i]) % f.size();

                // Get the meshpoint on my side
                label meshPointI = f[index];

                label localPointI = meshToLocalPoint
                (
                    meshToPatchPoint,
                    meshPointI
                );

                if (storeInfo(nbrInfo[i], localPointI))
                {
                    changedPoints.insert(localPointI);
                }
            }
        }
        else if (isA<cyclicPolyPatch>(pp))
        {
            // Handle cyclics: send lower half to upper half and vice versa.
            // Or since they both are in memory just do it point by point.

            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(pp);

            const labelList& meshPoints = pp.meshPoints();

            //const edgeList& connections = cycPatch.coupledPoints();
            const edgeList connections(coupledPoints(cycPatch));

            forAll(connections, i)
            {
                const edge& e = connections[i];

                label meshPointA = meshPoints[e[0]];
                label meshPointB = meshPoints[e[1]];

                label localA = meshToLocalPoint
                (
                    meshToPatchPoint,
                    meshPointA
                );
                label localB = meshToLocalPoint
                (
                    meshToPatchPoint,
                    meshPointB
                );


                // Do we have information on pointA?
                Map<label>::iterator procPointA =
                    meshToProcPoint_.find(localA);

                if (procPointA != meshToProcPoint_.end())
                {
                    // Store A info onto pointB
                    if (storeInfo(procPoints_[procPointA()], localB))
                    {
                        changedPoints.insert(localB);
                    }
                }

                // Same for info on pointB
                Map<label>::iterator procPointB =
                    meshToProcPoint_.find(localB);

                if (procPointB != meshToProcPoint_.end())
                {
                    // Store B info onto pointA
                    if (storeInfo(procPoints_[procPointB()], localA))
                    {
                        changedPoints.insert(localA);
                    }
                }
            }
        }
    }
}


// Remove entries which are handled by normal face-face communication. I.e.
// those points where the equivalence list is only me and my (face)neighbour
void Foam::globalPoints::remove
(
    const labelList& patchToMeshPoint,
    const Map<label>& directNeighbours
)
{
    // Save old ones.
    Map<label> oldMeshToProcPoint(meshToProcPoint_);
    meshToProcPoint_.clear();

    List<labelList> oldProcPoints;
    oldProcPoints.transfer(procPoints_);

    // Go through all equivalences
    forAllConstIter(Map<label>, oldMeshToProcPoint, iter)
    {
        label localPointI = iter.key();
        const labelList& pointInfo = oldProcPoints[iter()];

        if (pointInfo.size() == 2)
        {
            // I will be in this equivalence list.
            // Check whether my direct (=face) neighbour
            // is in it. This would be an ordinary connection and can be
            // handled by normal face-face connectivity.

            const label a = pointInfo[0];
            const label b = pointInfo[1];

            if
            (
                (
                    globalIndices_.isLocal(a)
                 && directNeighbours.found(globalIndices_.toLocal(a))
                )
             || (
                    globalIndices_.isLocal(b)
                 && directNeighbours.found(globalIndices_.toLocal(b))
                )
            )
            {
                // Normal faceNeighbours
                if (globalIndices_.isLocal(a))
                {
                    //Pout<< "Removing direct neighbour:"
                    //    << mesh_.points()[a[1]]
                    //    << endl;
                }
                else if (globalIndices_.isLocal(b))
                {
                    //Pout<< "Removing direct neighbour:"
                    //    << mesh_.points()[b[1]]
                    //    << endl;
                }
            }
            else
            {
                // This condition will be very rare: points are used by
                // two processors which are not face-face connected.
                // e.g.
                // +------+------+
                // | wall |  B   |
                // +------+------+
                // |   A  | wall |
                // +------+------+
                // Processor A and B share a point. Note that this only will
                // be found if the two domains are face connected at all
                // (not shown in the picture)

                meshToProcPoint_.insert(localPointI, procPoints_.size());
                procPoints_.append(pointInfo);
            }
        }
        else if (pointInfo.size() == 1)
        {
            // This happens for 'wedge' like cyclics where the two halves
            // come together in the same point so share the same meshPoint.
            // So this meshPoint will have info of size one only.
            if
            (
                !globalIndices_.isLocal(pointInfo[0])
             || !directNeighbours.found(globalIndices_.toLocal(pointInfo[0]))
            )
            {
                meshToProcPoint_.insert(localPointI, procPoints_.size());
                procPoints_.append(pointInfo);
            }
        }
        else
        {
            meshToProcPoint_.insert(localPointI, procPoints_.size());
            procPoints_.append(pointInfo);
        }
    }

    procPoints_.shrink();
}


// Compact indices
void Foam::globalPoints::compact()
{
    // TBD: find same procPoints entries. Or rather check if
    // in a procPoints there are points with indices < my index.
    // This will just merge identical entries so lower storage, but will
    // not affect anything else. Note: only relevant with cyclics.

    labelList oldToNew(procPoints_.size(), -1);
    labelList newToOld(meshToProcPoint_.size());

    label newIndex = 0;
    forAllIter(Map<label>, meshToProcPoint_, iter)
    {
        label oldIndex = iter();

        if (oldToNew[oldIndex] == -1)
        {
            iter() = newIndex;
            oldToNew[oldIndex] = newIndex;
            newToOld[newIndex] = oldIndex;
            newIndex++;
        }
    }
    List<labelList> oldProcPoints;
    oldProcPoints.transfer(procPoints_);

    procPoints_.setSize(meshToProcPoint_.size());
    forAll(procPoints_, i)
    {
        // Transfer
        procPoints_[i].transfer(oldProcPoints[newToOld[i]]);
    }
}


// Get (indices of) points for which I am master (= lowest numbered point on
// lowest numbered processor).
// (equivalence lists should be complete by now)
Foam::labelList Foam::globalPoints::getMasterPoints
(
    const labelList& patchToMeshPoint
) const
{
    labelList masterPoints(nPatchPoints_);
    label nMaster = 0;

    // Go through all equivalences and determine points where I am master.
    forAllConstIter(Map<label>, meshToProcPoint_, iter)
    {
        label localPointI = iter.key();
        const labelList& pointInfo = procPoints_[iter()];

        if (pointInfo.size() < 2)
        {
            // Points should have an equivalence list >= 2 since otherwise
            // they would be direct neighbours and have been removed in the
            // call to 'remove'.
            label meshPointI = localToMeshPoint(patchToMeshPoint, localPointI);

            FatalErrorIn("globalPoints::getMasterPoints(..)")
                << '[' << Pstream::myProcNo() << ']'
                << " MeshPoint:" << meshPointI
                << " coord:" << mesh_.points()[meshPointI]
                << " has no corresponding point on a neighbouring processor"
                << abort(FatalError);
        }
        else
        {
            // Check if lowest numbered processor and point is
            // on this processor. Since already sorted is first element.

            if
            (
                globalIndices_.isLocal(pointInfo[0])
             && globalIndices_.toLocal(pointInfo[0]) == localPointI
            )
            {
                // I am lowest numbered processor and point. Add to my list.
                masterPoints[nMaster++] = localPointI;
            }
        }
    }

    masterPoints.setSize(nMaster);

    return masterPoints;
}


// Send subset of lists
void Foam::globalPoints::sendSharedPoints
(
    PstreamBuffers& pBufs,
    const labelList& changedIndices
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (Pstream::parRun() && isA<processorPolyPatch>(pp))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pp);

            UOPstream toNeighbour(procPatch.neighbProcNo(), pBufs);

            if (debug)
            {
                Pout<< "Sending to " << procPatch.neighbProcNo()
                    << "  changed sharedPoints info:"
                    << changedIndices.size() << endl;
            }

            // Send over changed elements
            toNeighbour
                << UIndirectList<label>(sharedPointAddr_, changedIndices)()
                << UIndirectList<label>(sharedPointLabels_, changedIndices)();
        }
    }
}


// Receive shared point indices for all my shared points. Note that since
// there are only a few here we can build a reverse map using the point label
// instead of doing all this relative point indexing (patch face + index in
// face) as in send/receivePatchPoints
void Foam::globalPoints::receiveSharedPoints
(
    const Map<label>& meshToPatchPoint,
    PstreamBuffers& pBufs,
    labelList& changedIndices
)
{
    changedIndices.setSize(sharedPointAddr_.size());
    label nChanged = 0;

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Receive and set shared points
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (Pstream::parRun() && isA<processorPolyPatch>(pp))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pp);

            // Map from neighbouring mesh or patch point to sharedPoint)
            Map<label> nbrSharedPoints(sharedPointAddr_.size());

            {
                // Receive points on neighbour and sharedPoints and build
                // map from it. Note that we could have built the map on the
                // neighbour and sent it over.
                labelList nbrSharedPointAddr;
                labelList nbrSharedPointLabels;

                {
                    UIPstream fromNeighbour(procPatch.neighbProcNo(), pBufs);
                    fromNeighbour >> nbrSharedPointAddr >> nbrSharedPointLabels;
                }

                // Insert into to map
                forAll(nbrSharedPointLabels, i)
                {
                    nbrSharedPoints.insert
                    (
                        nbrSharedPointLabels[i], // mesh/patchpoint on neighbour
                        nbrSharedPointAddr[i]    // sharedPoint label
                    );
                }
            }


            // Merge into whatever information I hold.
            forAllConstIter(Map<label>, meshToProcPoint_, iter)
            {
                label localPointI = iter.key();
                label index = iter();

                if (sharedPointAddr_[index] == -1)
                {
                    // No shared point known yet for this point.
                    // See if was received from neighbour.
                    const labelList& knownInfo = procPoints_[index];

                    // Check through the whole equivalence list for any
                    // point from the neighbour.
                    forAll(knownInfo, j)
                    {
                        const label info = knownInfo[j];
                        label procI = globalIndices_.whichProcID(info);
                        label pointI = globalIndices_.toLocal(procI, info);

                        if
                        (
                            procI == procPatch.neighbProcNo()
                         && nbrSharedPoints.found(pointI)
                        )
                        {
                            // So this knownInfo contains the neighbour point
                            label sharedPointI = nbrSharedPoints[pointI];

                            sharedPointAddr_[index] = sharedPointI;
                            sharedPointLabels_[index] = localPointI;
                            changedIndices[nChanged++] = index;

                            break;
                        }
                    }
                }
            }
        }
        else if (isA<cyclicPolyPatch>(pp))
        {
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(pp);

            // Build map from mesh or patch point to sharedPoint
            Map<label> localToSharedPoint(sharedPointAddr_.size());
            forAll(sharedPointLabels_, i)
            {
                localToSharedPoint.insert
                (
                    sharedPointLabels_[i],
                    sharedPointAddr_[i]
                );
            }

            // Sync all info.
            //const edgeList& connections = cycPatch.coupledPoints();
            const edgeList connections(coupledPoints(cycPatch));

            forAll(connections, i)
            {
                const edge& e = connections[i];
                label meshPointA = pp.meshPoints()[e[0]];
                label meshPointB = pp.meshPoints()[e[1]];

                label localA = meshToLocalPoint
                (
                    meshToPatchPoint,
                    meshPointA
                );
                label localB = meshToLocalPoint
                (
                    meshToPatchPoint,
                    meshPointB
                );

                // Do we already have shared point for pointA?
                Map<label>::iterator fndA = localToSharedPoint.find(localA);
                Map<label>::iterator fndB = localToSharedPoint.find(localB);

                if (fndA != localToSharedPoint.end())
                {
                    if (fndB != localToSharedPoint.end())
                    {
                        if (fndA() != fndB())
                        {
                            FatalErrorIn
                            (
                                "globalPoints::receiveSharedPoints"
                                "(labelList&)"
                            )   << "On patch " << pp.name()
                                << " connected points " << meshPointA
                                << ' ' << mesh_.points()[meshPointA]
                                << " and " << meshPointB
                                << ' ' << mesh_.points()[meshPointB]
                                << " are mapped to different shared points: "
                                << fndA() << " and " << fndB()
                                << abort(FatalError);
                        }
                    }
                    else
                    {
                        // No shared point yet for B.
                        label sharedPointI = fndA();

                        // Store shared point for pointB
                        label index = meshToProcPoint_[localB];

                        sharedPointAddr_[index] = sharedPointI;
                        sharedPointLabels_[index] = localB;
                        changedIndices[nChanged++] = index;
                    }
                }
                else
                {
                    // No shared point yet for A.
                    if (fndB != localToSharedPoint.end())
                    {
                        label sharedPointI = fndB();

                        // Store shared point for pointA
                        label index = meshToProcPoint_[localA];

                        sharedPointAddr_[index] = sharedPointI;
                        sharedPointLabels_[index] = localA;
                        changedIndices[nChanged++] = index;
                    }
                }
            }
        }
    }

    changedIndices.setSize(nChanged);
}


Foam::edgeList Foam::globalPoints::coupledPoints(const cyclicPolyPatch& pp)
{
    // Look at cyclic patch as two halves, A and B.
    // Now all we know is that relative face index in halfA is same
    // as coupled face in halfB and also that the 0th vertex
    // corresponds.

    // From halfA point to halfB or -1.
    labelList coupledPoint(pp.nPoints(), -1);

    for (label patchFaceA = 0; patchFaceA < pp.size()/2; patchFaceA++)
    {
        const face& fA = pp.localFaces()[patchFaceA];

        forAll(fA, indexA)
        {
            label patchPointA = fA[indexA];

            if (coupledPoint[patchPointA] == -1)
            {
                const face& fB = pp.localFaces()[patchFaceA + pp.size()/2];

                label indexB = (fB.size() - indexA) % fB.size();

                coupledPoint[patchPointA] = fB[indexB];
            }
        }
    }

    edgeList connected(pp.nPoints());

    // Extract coupled points.
    label connectedI = 0;

    forAll(coupledPoint, i)
    {
        if (coupledPoint[i] != -1)
        {
            connected[connectedI++] = edge(i, coupledPoint[i]);
        }
    }

    connected.setSize(connectedI);

    return connected;
}


void Foam::globalPoints::calculateSharedPoints
(
    const Map<label>& meshToPatchPoint, // from mesh point to local numbering
    const labelList& patchToMeshPoint,  // from local numbering to mesh point
    const bool keepAllPoints
)
{
    if (debug)
    {
        Pout<< "globalPoints::globalPoints(const polyMesh&) : "
            << "doing processor to processor communication to get sharedPoints"
            << endl;
    }

    labelHashSet changedPoints(nPatchPoints_);

    // Initialize procPoints with my patch points. Keep track of points
    // inserted (in changedPoints)
    // There are two possible forms of this:
    // - initialize with all patch points (allPoints = true). This causes all
    //   patch points to be exchanged so a lot of information gets stored and
    //   transferred. This all gets filtered out later when removing the
    //   equivalence lists of size 2.
    // - initialize with boundary points of patches only (allPoints = false).
    //   This should work for all decompositions except extreme ones where a
    //   shared point is not on the boundary of any processor patches using it.
    //   This would happen if a domain was pinched such that two patches share
    //   a point or edge.
    initOwnPoints(meshToPatchPoint, true, changedPoints);

    // Do one exchange iteration to get neighbour points.
    {
        PstreamBuffers pBufs(Pstream::defaultCommsType);
        sendPatchPoints(meshToPatchPoint, pBufs, changedPoints);
        pBufs.finishedSends();
        receivePatchPoints(meshToPatchPoint, pBufs, changedPoints);
    }


    // Save neighbours reachable through face-face communication.
    Map<label> neighbourList;
    if (!keepAllPoints)
    {
        neighbourList = meshToProcPoint_;
    }

    // Exchange until nothing changes on all processors.
    bool changed = false;

    do
    {
        PstreamBuffers pBufs(Pstream::defaultCommsType);
        sendPatchPoints(meshToPatchPoint, pBufs, changedPoints);
        pBufs.finishedSends();
        receivePatchPoints(meshToPatchPoint, pBufs, changedPoints);

        changed = changedPoints.size() > 0;
        reduce(changed, orOp<bool>());

    } while (changed);


    // Remove direct neighbours from point equivalences.
    if (!keepAllPoints)
    {
        remove(patchToMeshPoint, neighbourList);
    }
    else
    {
        // Compact out unused elements
        compact();
    }

    procPoints_.shrink();

    // Sort procPoints in incremental order. This will make the master the
    // first element.
    forAllIter(Map<label>, meshToProcPoint_, iter)
    {
        sort(procPoints_[iter()]);
    }


    // Pout<< "Now connected points:" << endl;
    // forAllConstIter(Map<label>, meshToProcPoint_, iter)
    // {
    //     label localI = iter.key();
    //     const labelList& pointInfo = procPoints_[iter()];
    // 
    //     Pout<< "pointI:" << localI << " index:" << iter()
    //         << " coord:"
    //         << mesh_.points()[localToMeshPoint(patchToMeshPoint, localI)]
    //         << endl;
    // 
    //     forAll(pointInfo, i)
    //     {
    //         label procI = globalIndices_.whichProcID(pointInfo[i]);
    //         Pout<< "    connected to proc " << procI
    //             << " localpoint:"
    //             << globalIndices_.toLocal(procI, pointInfo[i]);
    // 
    //         if (globalIndices_.isLocal(pointInfo[i]))
    //         {
    //             label meshPointI = localToMeshPoint
    //             (
    //                 patchToMeshPoint,
    //                 globalIndices_.toLocal(pointInfo[i])
    //             );
    //             Pout<< " at:" <<  mesh_.points()[meshPointI];
    //         }
    // 
    //         Pout<< endl;
    //     }
    // }


    // We now have - in procPoints_ - a list of points which are shared between
    // multiple processors. These are the ones for which are sharedPoint
    // needs to be determined. This is done by having the lowest numbered
    // processor in the equivalence list 'ask' for a sharedPoint number
    // and then distribute it across processor patches to the non-master
    // processors. Note: below piece of coding is not very efficient. Uses
    // a Map where possibly it shouldn't

    // Initialize sharedPoint addressing. Is for every entry in procPoints_
    // the sharedPoint.
    sharedPointAddr_.setSize(meshToProcPoint_.size());
    sharedPointAddr_ = -1;
    sharedPointLabels_.setSize(meshToProcPoint_.size());
    sharedPointLabels_ = -1;


    // Get point labels of points for which I am master (lowest
    // numbered proc)
    labelList masterPoints(getMasterPoints(patchToMeshPoint));

    // Get global numbering for master points
    globalIndex globalMasterPoints(masterPoints.size());
    nGlobalPoints_ = globalMasterPoints.size();

    forAll(masterPoints, i)
    {
        label localPointI = masterPoints[i];
        label index = meshToProcPoint_[localPointI];

        sharedPointLabels_[index] = localPointI;
        sharedPointAddr_[index] = globalMasterPoints.toGlobal(i);
    }


    // Now we have a sharedPointLabel for some of the entries in procPoints.
    // Send this information to neighbours. Receive their information.
    // Loop until nothing changes.

    // Initial subset to send is points for which I have sharedPoints
    labelList changedIndices(sharedPointAddr_.size());
    label nChanged = 0;

    forAll(sharedPointAddr_, i)
    {
        if (sharedPointAddr_[i] != -1)
        {
            changedIndices[nChanged++] = i;
        }
    }
    changedIndices.setSize(nChanged);

    changed = false;

    do
    {
        if (debug)
        {
            Pout<< "Determined " << changedIndices.size() << " shared points."
                << " Exchanging them" << endl;
        }
        PstreamBuffers pBufs(Pstream::defaultCommsType);
        sendSharedPoints(pBufs, changedIndices);
        pBufs.finishedSends();
        receiveSharedPoints(meshToPatchPoint, pBufs, changedIndices);

        changed = changedIndices.size() > 0;
        reduce(changed, orOp<bool>());

    } while (changed);


    forAll(sharedPointLabels_, i)
    {
        if (sharedPointLabels_[i] == -1)
        {
            FatalErrorIn("globalPoints::globalPoints(const polyMesh&)")
                << "Problem: shared point on processor " << Pstream::myProcNo()
                << " not set at index " << sharedPointLabels_[i] << endl
                << "This might mean the individual processor domains are not"
                << " connected and the overall domain consists of multiple"
                << " regions. You can check this with checkMesh"
                << abort(FatalError);
        }
    }

    if (debug)
    {
        Pout<< "globalPoints::globalPoints(const polyMesh&) : "
            << "Finished global points" << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::globalPoints::globalPoints
(
    const polyMesh& mesh,
    const bool keepAllPoints
)
:
    mesh_(mesh),
    globalIndices_(mesh_.nPoints()),
    nPatchPoints_(countPatchPoints(mesh.boundaryMesh())),
    procPoints_(nPatchPoints_),
    meshToProcPoint_(nPatchPoints_),
    sharedPointAddr_(0),
    sharedPointLabels_(0),
    nGlobalPoints_(0)
{
    // Empty patch maps to signal storing mesh point labels
    Map<label> meshToPatchPoint(0);
    labelList patchToMeshPoint(0);

    calculateSharedPoints(meshToPatchPoint, patchToMeshPoint, keepAllPoints);
}


// Construct from mesh and patch of coupled faces
Foam::globalPoints::globalPoints
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& coupledPatch,
    const bool keepAllPoints
)
:
    mesh_(mesh),
    globalIndices_(coupledPatch.nPoints()),
    nPatchPoints_(coupledPatch.nPoints()),
    procPoints_(nPatchPoints_),
    meshToProcPoint_(nPatchPoints_),
    sharedPointAddr_(0),
    sharedPointLabels_(0),
    nGlobalPoints_(0)
{
    calculateSharedPoints
    (
        coupledPatch.meshPointMap(),
        coupledPatch.meshPoints(),
        keepAllPoints
    );
}


// ************************************************************************* //
