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

const Foam::label Foam::globalPoints::fromCollocated = labelMax/2;

const Foam::scalar Foam::globalPoints::mergeDist = ROOTVSMALL;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Routines to handle global indices
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bool Foam::globalPoints::noTransform(const tensor& tt, const scalar mergeDist)
{
    return
        (mag(tt.xx()-1) < mergeDist)
     && (mag(tt.yy()-1) < mergeDist)
     && (mag(tt.zz()-1) < mergeDist)
     && (mag(tt.xy()) < mergeDist)
     && (mag(tt.xz()) < mergeDist)
     && (mag(tt.yx()) < mergeDist)
     && (mag(tt.yz()) < mergeDist)
     && (mag(tt.zx()) < mergeDist)
     && (mag(tt.zy()) < mergeDist);
}


// Calculates per face whether couple is collocated.
Foam::PackedBoolList Foam::globalPoints::collocatedFaces
(
    const coupledPolyPatch& pp,
    const scalar mergeDist
)
{
    // Initialise to false
    PackedBoolList collocated(pp.size());

    const vectorField& separation = pp.separation();
    const tensorField& forwardT = pp.forwardT();

    if (forwardT.size() == 0)
    {
        // Parallel.
        if (separation.size() == 0)
        {
            collocated = 1u;
        }
        else if (separation.size() == 1)
        {
            // Fully separate. Do not synchronise.
        }
        else
        {
            // Per face separation.
            forAll(pp, faceI)
            {
                if (mag(separation[faceI]) < mergeDist)
                {
                    collocated[faceI] = 1u;
                }
            }
        }
    }
    else if (forwardT.size() == 1)
    {
        // Fully transformed.
    }
    else
    {
        // Per face transformation.
        forAll(pp, faceI)
        {
            if (noTransform(forwardT[faceI], mergeDist))
            {
                collocated[faceI] = 1u;
            }
        }
    }
    return collocated;
}


Foam::PackedBoolList Foam::globalPoints::collocatedPoints
(
    const coupledPolyPatch& pp,
    const scalar mergeDist
)
{
    // Initialise to false
    PackedBoolList collocated(pp.nPoints());

    const vectorField& separation = pp.separation();
    const tensorField& forwardT = pp.forwardT();

    if (forwardT.size() == 0)
    {
        // Parallel.
        if (separation.size() == 0)
        {
            collocated = 1u;
        }
        else if (separation.size() == 1)
        {
            // Fully separate.
        }
        else
        {
            // Per face separation.
            for (label pointI = 0; pointI < pp.nPoints(); pointI++)
            {
                label faceI = pp.pointFaces()[pointI][0];

                if (mag(separation[faceI]) < mergeDist)
                {
                    collocated[pointI] = 1u;
                }
            }
        }
    }
    else if (forwardT.size() == 1)
    {
        // Fully transformed.
    }
    else
    {
        // Per face transformation.
        for (label pointI = 0; pointI < pp.nPoints(); pointI++)
        {
            label faceI = pp.pointFaces()[pointI][0];

            if (noTransform(forwardT[faceI], mergeDist))
            {
                collocated[pointI] = 1u;
            }
        }
    }
    return collocated;
}


Foam::label Foam::globalPoints::toGlobal
(
    const label localPointI,
    const bool isCollocated
) const
{
    label globalPointI = globalIndices_.toGlobal(localPointI);

    if (isCollocated)
    {
        return globalPointI + fromCollocated;
    }
    else
    {
        return globalPointI;
    }
}


bool Foam::globalPoints::isCollocated(const label globalI) const
{
    return globalI >= fromCollocated;
}


Foam::label Foam::globalPoints::removeCollocated(const label globalI) const
{
    if (globalI >= fromCollocated)
    {
        return globalI - fromCollocated;
    }
    else
    {
        return globalI;
    }
}


bool Foam::globalPoints::isLocal(const label globalI) const
{
    return globalIndices_.isLocal(removeCollocated(globalI));
}


Foam::label Foam::globalPoints::whichProcID(const label globalI) const
{
    return globalIndices_.whichProcID(removeCollocated(globalI));
}


Foam::label Foam::globalPoints::toLocal
(
    const label procI,
    const label globalI
) const
{
    return globalIndices_.toLocal(procI, removeCollocated(globalI));
}


Foam::label Foam::globalPoints::toLocal(const label globalI) const
{
    return toLocal(Pstream::myProcNo(), globalI);
}


// Collect all topological information about a point on a patch.

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
    const label localPointI,
    const bool isCollocated
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
        labelList knownInfo(1, toGlobal(localPointI, isCollocated));

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


void Foam::globalPoints::printProcPoints
(
    const labelList& patchToMeshPoint,
    const labelList& pointInfo,
    Ostream& os
) const
{
    forAll(pointInfo, i)
    {
        label globalI = pointInfo[i];

        label procI = whichProcID(globalI);
        os  << "    connected to proc " << procI;

        if (isCollocated(globalI))
        {
            os  << " collocated localpoint:";
        }
        else
        {
            os  << " separated localpoint:";
        }

        os  << toLocal(procI, globalI);

        if (isLocal(globalI))
        {
            label meshPointI = localToMeshPoint
            (
                patchToMeshPoint,
                toLocal(globalI)
            );
            os  << " at:" <<  mesh_.points()[meshPointI];
        }

        os  << endl;
    }
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
            // Find points with transforms

            PackedBoolList isCollocatedPoint
            (
                collocatedPoints
                (
                    refCast<const coupledPolyPatch>(pp),
                    mergeDist
                )
            );


            const labelList& meshPoints = pp.meshPoints();

            if (allPoints)
            {
                // All points on patch
                forAll(meshPoints, patchPointI)
                {
                    label meshPointI = meshPoints[patchPointI];
                    label localPointI = meshToLocalPoint
                    (
                        meshToPatchPoint,
                        meshPointI
                    );
                    labelList knownInfo
                    (
                        1,
                        toGlobal(localPointI, isCollocatedPoint[patchPointI])
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
                        toGlobal
                        (
                            localPointI,
                            isCollocatedPoint[boundaryPoints[i]]
                        )
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
    const bool mergeSeparated,
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
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pp);

            PackedBoolList isCollocatedPoint
            (
                collocatedPoints
                (
                    procPatch,
                    mergeDist
                )
            );


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
                if (mergeSeparated || isCollocatedPoint[patchPointI])
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

                        // Add my information about localPointI to the
                        // send buffers
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
            }

            // Send to neighbour
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


// Receive all my neighbours' information and merge with mine.
// After finishing will have updated
// - procPoints_ : all neighbour information merged in.
// - meshToProcPoint_
// - changedPoints: all points for which something changed.
void Foam::globalPoints::receivePatchPoints
(
    const bool mergeSeparated,
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

            PackedBoolList isCollocatedPoint
            (
                collocatedPoints
                (
                    procPatch,
                    mergeDist
                )
            );

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

                if
                (
                    storeInfo
                    (
                        nbrInfo[i],
                        localPointI,
                        isCollocatedPoint[pp.meshPointMap()[meshPointI]]
                    )
                )
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

            PackedBoolList isCollocatedPoint
            (
                collocatedPoints
                (
                    cycPatch,
                    mergeDist
                )
            );

            const labelList& meshPoints = pp.meshPoints();

            //const edgeList& connections = cycPatch.coupledPoints();
            const edgeList connections(coupledPoints(cycPatch));

            forAll(connections, i)
            {
                const edge& e = connections[i];

                if (mergeSeparated || isCollocatedPoint[e[0]])
                {
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
                        if
                        (
                            storeInfo
                            (
                                procPoints_[procPointA()],
                                localB,
                                isCollocatedPoint[e[1]]
                            )
                        )
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
                        if
                        (
                            storeInfo
                            (
                                procPoints_[procPointB()],
                                localA,
                                isCollocatedPoint[e[0]]
                            )
                        )
                        {
                            changedPoints.insert(localA);
                        }
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
                    isLocal(a)
                 && directNeighbours.found(toLocal(a))
                )
             || (
                    isLocal(b)
                 && directNeighbours.found(toLocal(b))
                )
            )
            {
                // Normal faceNeighbours
                if (isLocal(a))
                {
                    //Pout<< "Removing direct neighbour:"
                    //    << mesh_.points()[a[1]]
                    //    << endl;
                }
                else if (isLocal(b))
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
                !isLocal(pointInfo[0])
             || !directNeighbours.found(toLocal(pointInfo[0]))
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
void Foam::globalPoints::compact(const labelList& patchToMeshPoint)
{
    labelList oldToNew(procPoints_.size(), -1);
    labelList newToOld(meshToProcPoint_.size());

    label newIndex = 0;
    forAllConstIter(Map<label>, meshToProcPoint_, iter)
    {
        label oldIndex = iter();

        if (oldToNew[oldIndex] == -1)
        {
            // Check if one of the procPoints already is merged.
            const labelList& pointInfo = procPoints_[oldIndex];

            if (pointInfo.size() >= 2)
            {
                // Merge out single entries

                label minIndex = labelMax;

                forAll(pointInfo, i)
                {
                    if (isLocal(pointInfo[i]))
                    {
                        label localI = toLocal(pointInfo[i]);

                        // Is localPoint itself already merged?
                        label index = meshToProcPoint_[localI];

                        //Pout<< "        found point:" << localI
                        //    << " with current index " << index << endl;

                        if (oldToNew[index] != -1)
                        {
                            minIndex = min(minIndex, oldToNew[index]);
                        }
                    }
                }

                if (minIndex < labelMax && minIndex != oldIndex)
                {
                    // Make my index point to minIndex
                    oldToNew[oldIndex] = minIndex;
                }
                else
                {
                    // Nothing compacted yet. Allocate new index.
                    oldToNew[oldIndex] = newIndex;
                    newToOld[newIndex] = oldIndex;
                    newIndex++;
                }
            }
            else
            {
                //Pout<< "Removed singlepoint pointI:" << iter.key()
                //    << " currentindex:" << iter()
                //    << endl;
            }
        }
    }

    // Redo indices
    Map<label> oldMeshToProcPoint(meshToProcPoint_.xfer());
    forAllConstIter(Map<label>, oldMeshToProcPoint, iter)
    {
        label newIndex = oldToNew[iter()];
        if (newIndex != -1)
        {
            meshToProcPoint_.insert(iter.key(), newIndex);
        }
    }

    List<labelList> oldProcPoints;
    oldProcPoints.transfer(procPoints_);

    newToOld.setSize(newIndex);

    procPoints_.setSize(newIndex);

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
    //labelList masterPoints(nPatchPoints_);
    //label nMaster = 0;

    labelHashSet masterPointSet(nPatchPoints_);


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
                isLocal(pointInfo[0])
             && toLocal(pointInfo[0]) == localPointI
            )
            {
                // I am lowest numbered processor and point. Add to my list.
                //masterPoints[nMaster++] = localPointI;
                masterPointSet.insert(localPointI);
            }
        }
    }

    return masterPointSet.toc();
}


// Send subset of lists.
// Note: might not be used if separated
void Foam::globalPoints::sendSharedPoints
(
    const bool mergeSeparated,
    PstreamBuffers& pBufs,
    const DynamicList<label>& changedIndices
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


// Check if any connected local points already have a sharedPoint. This handles
// locally connected points.
void Foam::globalPoints::extendSharedPoints
(
    const Map<label>& meshToShared,
    DynamicList<label>& changedIndices
)
{
    forAllConstIter(Map<label>, meshToProcPoint_, iter)
    {
        label localI = iter.key();
        label sharedI = meshToShared[localI];

        if (sharedPointLabels_[sharedI] == -1)
        {
            // No sharedpoint yet for me. Check if any of my connected
            // points have.

            const labelList& pointInfo = procPoints_[iter()];
            forAll(pointInfo, i)
            {
                if (isLocal(pointInfo[i]))
                {
                    label nbrLocalI = toLocal(pointInfo[i]);
                    label nbrSharedI = meshToShared[nbrLocalI];

                    if (sharedPointAddr_[nbrSharedI] != -1)
                    {
                        // I do have a sharedpoint for nbrSharedI. Reuse it.
                        sharedPointLabels_[sharedI] = localI;
                        sharedPointAddr_[sharedI] =
                            sharedPointAddr_[nbrSharedI];
                        changedIndices.append(sharedI);
                        break;
                    }
                }
            }
        }
    }
}


// Receive shared point indices for all my shared points. Note that since
// there are only a few here we can build a reverse map using the point label
// instead of doing all this relative point indexing (patch face + index in
// face) as in send/receivePatchPoints
void Foam::globalPoints::receiveSharedPoints
(
    const bool mergeSeparated,
    const Map<label>& meshToPatchPoint,
    const Map<label>& meshToShared,
    PstreamBuffers& pBufs,
    DynamicList<label>& changedIndices
)
{
    changedIndices.clear();

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
                label sharedI = meshToShared[localPointI];

                if (sharedPointAddr_[sharedI] == -1)
                {
                    // No shared point known yet for this point.
                    // See if was received from neighbour.
                    const labelList& knownInfo = procPoints_[iter()];

                    // Check through the whole equivalence list for any
                    // point from the neighbour.
                    forAll(knownInfo, j)
                    {
                        const label info = knownInfo[j];
                        label procI = whichProcID(info);
                        label pointI = toLocal(procI, info);

                        if
                        (
                            procI == procPatch.neighbProcNo()
                         && nbrSharedPoints.found(pointI)
                        )
                        {
                            // So this knownInfo contains the neighbour point
                            label sharedPointI = nbrSharedPoints[pointI];

                            sharedPointAddr_[sharedI] = sharedPointI;
                            sharedPointLabels_[sharedI] = localPointI;
                            changedIndices.append(sharedI);

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

            PackedBoolList isCollocatedPoint
            (
                collocatedPoints
                (
                    cycPatch,
                    mergeDist
                )
            );

            // Build map from mesh or patch point to sharedPoint.
            Map<label> localToSharedPoint(sharedPointAddr_.size());
            forAll(sharedPointLabels_, i)
            {
                if (sharedPointLabels_[i] != -1)
                {
                    localToSharedPoint.insert
                    (
                        sharedPointLabels_[i],
                        sharedPointAddr_[i]
                    );
                }
            }

            // Sync all info.
            //const edgeList& connections = cycPatch.coupledPoints();
            const edgeList connections(coupledPoints(cycPatch));

            forAll(connections, i)
            {
                const edge& e = connections[i];

                if (mergeSeparated || isCollocatedPoint[e[0]])
                {
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
                                    "globalPoints::receiveSharedPoints(..)"
                                )   << "On patch " << pp.name()
                                    << " connected points " << meshPointA
                                    << ' ' << mesh_.points()[meshPointA]
                                    << " and " << meshPointB
                                    << ' ' << mesh_.points()[meshPointB]
                                    << " are mapped to different shared"
                                    << " points: "
                                    << fndA() << " and " << fndB()
                                    << abort(FatalError);
                            }
                        }
                        else
                        {
                            // No shared point yet for B.
                            label sharedPointI = fndA();

                            // Store shared point for pointB
                            label sharedI = meshToShared[localB];

                            sharedPointAddr_[sharedI] = sharedPointI;
                            sharedPointLabels_[sharedI] = localB;
                            changedIndices.append(sharedI);
                        }
                    }
                    else
                    {
                        // No shared point yet for A.
                        if (fndB != localToSharedPoint.end())
                        {
                            label sharedPointI = fndB();

                            // Store shared point for pointA
                            label sharedI = meshToShared[localA];

                            sharedPointAddr_[sharedI] = sharedPointI;
                            sharedPointLabels_[sharedI] = localA;
                            changedIndices.append(sharedI);
                        }
                    }
                }
            }
        }
    }

    // Extend shared points from local connections.
    extendSharedPoints(meshToShared, changedIndices);
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
    const bool keepAllPoints,
    const bool mergeSeparated
)
{
    if (debug)
    {
        Pout<< "globalPoints::calculateSharedPoints(..) : "
            << "doing processor to processor communication to get sharedPoints"
            << endl;
    }

    if (globalIndices_.size() >= fromCollocated)
    {
        FatalErrorIn("globalPoints::calculateSharedPoints(..)")
            << "Integer overflow. Total number of points "
            << globalIndices_.size()
            << " is too large to be represented." << exit(FatalError);
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
        sendPatchPoints
        (
            mergeSeparated,
            meshToPatchPoint,
            pBufs,
            changedPoints
        );
        pBufs.finishedSends();
        receivePatchPoints
        (
            mergeSeparated,
            meshToPatchPoint,
            pBufs,
            changedPoints
        );
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
        sendPatchPoints
        (
            mergeSeparated,
            meshToPatchPoint,
            pBufs,
            changedPoints
        );
        pBufs.finishedSends();
        receivePatchPoints
        (
            mergeSeparated,
            meshToPatchPoint,
            pBufs,
            changedPoints
        );

        changed = changedPoints.size() > 0;
        reduce(changed, orOp<bool>());

    } while (changed);


    // Remove direct neighbours from point equivalences.
    if (!keepAllPoints)
    {
        //Pout<< "**Removing direct neighbours:" << endl;
        //forAllConstIter(Map<label>, neighbourList, iter)
        //{
        //    label localI = iter.key();
        //    Pout<< "pointI:" << localI << " index:" << iter()
        //        << " coord:"
        //        << mesh_.points()[localToMeshPoint(patchToMeshPoint, localI)]
        //        << endl;
        //}
        remove(patchToMeshPoint, neighbourList);
    }


    // Compact out unused elements
    compact(patchToMeshPoint);



    procPoints_.shrink();

    // Sort procPoints in incremental order. This will make the master the
    // first element.
    forAllIter(Map<label>, meshToProcPoint_, iter)
    {
        sort(procPoints_[iter()]);
    }
    // We can now remove the collocated bit
    forAll(procPoints_, index)
    {
        labelList& pointInfo = procPoints_[index];

        forAll(pointInfo, i)
        {
            pointInfo[i] = removeCollocated(pointInfo[i]);
        }
    }


//Pout<< "Now connected points:" << endl;
//forAllConstIter(Map<label>, meshToProcPoint_, iter)
//{
//    label localI = iter.key();
//    const labelList& pointInfo = procPoints_[iter()];
//
//    Pout<< "pointI:" << localI << " index:" << iter()
//        << " coord:"
//        << mesh_.points()[localToMeshPoint(patchToMeshPoint, localI)]
//        << endl;
//
//    printProcPoints(patchToMeshPoint, pointInfo, Pout);
//}


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

    // Number meshToProcPoint
    label sharedI = 0;
    Map<label> meshToShared(meshToProcPoint_.size());
    forAllConstIter(Map<label>, meshToProcPoint_, iter)
    {
        meshToShared.insert(iter.key(), sharedI++);
    }

    // Assign the entries for the masters
    forAll(masterPoints, i)
    {
        label localPointI = masterPoints[i];
        label sharedI = meshToShared[localPointI];

        sharedPointLabels_[sharedI] = localPointI;
        sharedPointAddr_[sharedI] = globalMasterPoints.toGlobal(i);
    }


    // Now we have a sharedPointLabel for some of the entries in procPoints.
    // Send this information to neighbours. Receive their information.
    // Loop until nothing changes.

    // Initial subset to send is points for which I have sharedPoints
    DynamicList<label> changedIndices(sharedPointAddr_.size());
    forAll(sharedPointAddr_, i)
    {
        if (sharedPointAddr_[i] != -1)
        {
            changedIndices.append(i);
        }
    }

    // Assign local slaves from local master
    extendSharedPoints(meshToShared, changedIndices);


    changed = false;

    do
    {
        if (debug)
        {
            Pout<< "Determined " << changedIndices.size() << " shared points."
                << " Exchanging them" << endl;
        }
        PstreamBuffers pBufs(Pstream::defaultCommsType);
        sendSharedPoints(mergeSeparated, pBufs, changedIndices);
        pBufs.finishedSends();
        receiveSharedPoints
        (
            mergeSeparated,
            meshToPatchPoint,
            meshToShared,
            pBufs,
            changedIndices
        );

        changed = changedIndices.size() > 0;
        reduce(changed, orOp<bool>());

    } while (changed);


//    forAll(sharedPointLabels_, sharedI)
//    {
//        label localI = sharedPointLabels_[sharedI];
//
//        Pout<< "point:" << localI
//            << " coord:"
//            << mesh_.points()[localToMeshPoint(patchToMeshPoint, localI)]
//            << " ON TO GLOBAL:" << sharedPointAddr_[sharedI]
//            << endl;
//    }


    forAll(sharedPointLabels_, i)
    {
        if (sharedPointLabels_[i] == -1)
        {
            FatalErrorIn("globalPoints::calculateSharedPoints(..)")
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
        Pout<< "globalPoints::calculateSharedPoints(..) : "
            << "Finished global points" << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::globalPoints::globalPoints
(
    const polyMesh& mesh,
    const bool keepAllPoints,
    const bool mergeSeparated
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

    calculateSharedPoints
    (
        meshToPatchPoint,
        patchToMeshPoint,
        keepAllPoints,
        mergeSeparated
    );
}


// Construct from mesh and patch of coupled faces
Foam::globalPoints::globalPoints
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& coupledPatch,
    const bool keepAllPoints,
    const bool mergeSeparated
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
        keepAllPoints,
        mergeSeparated
    );
}


// ************************************************************************* //
