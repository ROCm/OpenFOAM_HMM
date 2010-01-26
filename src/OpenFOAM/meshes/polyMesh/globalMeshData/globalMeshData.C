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

#include "globalMeshData.H"
#include "Time.H"
#include "Pstream.H"
#include "PstreamCombineReduceOps.H"
#include "processorPolyPatch.H"
#include "demandDrivenData.H"
#include "globalPoints.H"
//#include "geomGlobalPoints.H"
#include "polyMesh.H"
#include "mapDistribute.H"
#include "labelIOList.H"
#include "PackedList.H"
#include "mergePoints.H"
#include "matchPoints.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::globalMeshData, 0);

// Geometric matching tolerance. Factor of mesh bounding box.
const Foam::scalar Foam::globalMeshData::matchTol_ = 1E-8;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Collect processor patch addressing.
void Foam::globalMeshData::initProcAddr()
{
    processorPatchIndices_.setSize(mesh_.boundaryMesh().size());
    processorPatchIndices_ = -1;

    processorPatchNeighbours_.setSize(mesh_.boundaryMesh().size());
    processorPatchNeighbours_ = -1;

    // Construct processor patch indexing. processorPatchNeighbours_ only
    // set if running in parallel!
    processorPatches_.setSize(mesh_.boundaryMesh().size());

    label nNeighbours = 0;

    forAll (mesh_.boundaryMesh(), patchi)
    {
        if (isA<processorPolyPatch>(mesh_.boundaryMesh()[patchi]))
        {
            processorPatches_[nNeighbours] = patchi;
            processorPatchIndices_[patchi] = nNeighbours++;
        }
    }
    processorPatches_.setSize(nNeighbours);


    if (Pstream::parRun())
    {
        // Send indices of my processor patches to my neighbours
        forAll (processorPatches_, i)
        {
            label patchi = processorPatches_[i];

            OPstream toNeighbour
            (
                Pstream::blocking,
                refCast<const processorPolyPatch>
                (
                    mesh_.boundaryMesh()[patchi]
                ).neighbProcNo()
            );

            toNeighbour << processorPatchIndices_[patchi];
        }

        forAll(processorPatches_, i)
        {
            label patchi = processorPatches_[i];

            IPstream fromNeighbour
            (
                Pstream::blocking,
                refCast<const processorPolyPatch>
                (
                    mesh_.boundaryMesh()[patchi]
                ).neighbProcNo()
            );

            fromNeighbour >> processorPatchNeighbours_[patchi];
        }
    }
}


// Given information about locally used edges allocate global shared edges.
void Foam::globalMeshData::countSharedEdges
(
    const EdgeMap<labelList>& procSharedEdges,
    EdgeMap<label>& globalShared,
    label& sharedEdgeI
)
{
    // Count occurrences of procSharedEdges in global shared edges table.
    forAllConstIter(EdgeMap<labelList>, procSharedEdges, iter)
    {
        const edge& e = iter.key();

        EdgeMap<label>::iterator globalFnd = globalShared.find(e);

        if (globalFnd == globalShared.end())
        {
            // First time occurrence of this edge. Check how many we are adding.
            if (iter().size() == 1)
            {
                // Only one edge. Mark with special value.
                globalShared.insert(e, -1);
            }
            else
            {
                // Edge used more than once (even by local shared edges alone)
                // so allocate proper shared edge label.
                globalShared.insert(e, sharedEdgeI++);
            }
        }
        else
        {
            if (globalFnd() == -1)
            {
                // Second time occurence of this edge. Assign proper
                // edge label.
                globalFnd() = sharedEdgeI++;
            }
        }
    }
}


// Shared edges are shared between multiple processors. By their nature both
// of their endpoints are shared points. (but not all edges using two shared
// points are shared edges! There might e.g. be an edge between two unrelated
// clusters of shared points)
void Foam::globalMeshData::calcSharedEdges() const
{
    if (nGlobalEdges_ != -1 || sharedEdgeLabelsPtr_ || sharedEdgeAddrPtr_)
    {
        FatalErrorIn("globalMeshData::calcSharedEdges()")
            << "Shared edge addressing already done" << abort(FatalError);
    }


    const labelList& sharedPtAddr = sharedPointAddr();
    const labelList& sharedPtLabels = sharedPointLabels();

    // Since don't want to construct pointEdges for whole mesh create
    // Map for all shared points.
    Map<label> meshToShared(2*sharedPtLabels.size());
    forAll(sharedPtLabels, i)
    {
        meshToShared.insert(sharedPtLabels[i], i);
    }

    // Find edges using shared points. Store correspondence to local edge
    // numbering. Note that multiple local edges can have the same shared
    // points! (for cyclics or separated processor patches)
    EdgeMap<labelList> localShared(2*sharedPtAddr.size());

    const edgeList& edges = mesh_.edges();

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        Map<label>::const_iterator e0Fnd = meshToShared.find(e[0]);

        if (e0Fnd != meshToShared.end())
        {
            Map<label>::const_iterator e1Fnd = meshToShared.find(e[1]);

            if (e1Fnd != meshToShared.end())
            {
                // Found edge which uses shared points. Probably shared.

                // Construct the edge in shared points (or rather global indices
                // of the shared points)
                edge sharedEdge
                (
                    sharedPtAddr[e0Fnd()],
                    sharedPtAddr[e1Fnd()]
                );

                EdgeMap<labelList>::iterator iter =
                    localShared.find(sharedEdge);

                if (iter == localShared.end())
                {
                    // First occurrence of this point combination. Store.
                    localShared.insert(sharedEdge, labelList(1, edgeI));
                }
                else
                {
                    // Add this edge to list of edge labels.
                    labelList& edgeLabels = iter();

                    label sz = edgeLabels.size();
                    edgeLabels.setSize(sz+1);
                    edgeLabels[sz] = edgeI;
                }
            }
        }
    }


    // Now we have a table on every processors which gives its edges which use
    // shared points. Send this all to the master and have it allocate
    // global edge numbers for it. But only allocate a global edge number for
    // edge if it is used more than once!
    // Note that we are now sending the whole localShared to the master whereas
    // we only need the local count (i.e. the number of times a global edge is
    // used). But then this only gets done once so not too bothered about the
    // extra global communication.

    EdgeMap<label> globalShared(nGlobalPoints());

    if (Pstream::master())
    {
        label sharedEdgeI = 0;

        // Merge my shared edges into the global list
        if (debug)
        {
            Pout<< "globalMeshData::calcSharedEdges : Merging in from proc0 : "
                << localShared.size() << endl;
        }
        countSharedEdges(localShared, globalShared, sharedEdgeI);

        // Receive data from slaves and insert
        if (Pstream::parRun())
        {
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                // Receive the edges using shared points from the slave.
                IPstream fromSlave(Pstream::blocking, slave);
                EdgeMap<labelList> procSharedEdges(fromSlave);

                if (debug)
                {
                    Pout<< "globalMeshData::calcSharedEdges : "
                        << "Merging in from proc"
                        << Foam::name(slave) << " : " << procSharedEdges.size()
                        << endl;
                }
                countSharedEdges(procSharedEdges, globalShared, sharedEdgeI);
            }
        }

        // Now our globalShared should have some edges with -1 as edge label
        // These were only used once so are not proper shared edges.
        // Remove them.
        {
            EdgeMap<label> oldSharedEdges(globalShared);

            globalShared.clear();

            forAllConstIter(EdgeMap<label>, oldSharedEdges, iter)
            {
                if (iter() != -1)
                {
                    globalShared.insert(iter.key(), iter());
                }
            }
            if (debug)
            {
                Pout<< "globalMeshData::calcSharedEdges : Filtered "
                    << oldSharedEdges.size()
                    << " down to " << globalShared.size() << endl;
            }
        }


        // Send back to slaves.
        if (Pstream::parRun())
        {
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                // Receive the edges using shared points from the slave.
                OPstream toSlave(Pstream::blocking, slave);
                toSlave << globalShared;
            }
        }
    }
    else
    {
        // Send local edges to master
        {
            OPstream toMaster(Pstream::blocking, Pstream::masterNo());

            toMaster << localShared;
        }
        // Receive merged edges from master.
        {
            IPstream fromMaster(Pstream::blocking, Pstream::masterNo());

            fromMaster >> globalShared;
        }
    }

    // Now use the global shared edges list (globalShared) to classify my local
    // ones (localShared)

    nGlobalEdges_ = globalShared.size();

    DynamicList<label> dynSharedEdgeLabels(globalShared.size());
    DynamicList<label> dynSharedEdgeAddr(globalShared.size());

    forAllConstIter(EdgeMap<labelList>, localShared, iter)
    {
        const edge& e = iter.key();

        EdgeMap<label>::const_iterator edgeFnd = globalShared.find(e);

        if (edgeFnd != globalShared.end())
        {
            // My local edge is indeed a shared one. Go through all local edge
            // labels with this point combination.
            const labelList& edgeLabels = iter();

            forAll(edgeLabels, i)
            {
                // Store label of local mesh edge
                dynSharedEdgeLabels.append(edgeLabels[i]);

                // Store label of shared edge
                dynSharedEdgeAddr.append(edgeFnd());
            }
        }
    }

    sharedEdgeLabelsPtr_ = new labelList();
    labelList& sharedEdgeLabels = *sharedEdgeLabelsPtr_;
    sharedEdgeLabels.transfer(dynSharedEdgeLabels);

    sharedEdgeAddrPtr_ = new labelList();
    labelList& sharedEdgeAddr = *sharedEdgeAddrPtr_;
    sharedEdgeAddr.transfer(dynSharedEdgeAddr);

    if (debug)
    {
        Pout<< "globalMeshData : nGlobalEdges_:" << nGlobalEdges_ << nl
            << "globalMeshData : sharedEdgeLabels:" << sharedEdgeLabels.size()
            << nl
            << "globalMeshData : sharedEdgeAddr:" << sharedEdgeAddr.size()
            << endl;
    }
}


// Helper function to count coincident faces. This part used to be
// in updateMesh but I've had to move it to a separate function
// because of aliasing optimisation errors in icc9.1 on the
// Itanium.
Foam::label Foam::globalMeshData::countCoincidentFaces
(
    const scalar tolDim,
    const vectorField& separationDist
)
{
    label nCoincident = 0;

    forAll(separationDist, faceI)
    {
        if (mag(separationDist[faceI]) < tolDim)
        {
            // Faces coincide
            nCoincident++;
        }
    }
    return nCoincident;
}


void Foam::globalMeshData::calcGlobalPointSlaves() const
{
    if (debug)
    {
        Pout<< "globalMeshData::calcGlobalPointSlaves() :"
            << " calculating coupled master to slave point addressing."
            << endl;
    }

    // Calculate connected points for master points
    globalPoints globalData(mesh_, coupledPatch(), true);

    const Map<label>& meshToProcPoint = globalData.meshToProcPoint();

    // Create global numbering for coupled points
    globalPointNumberingPtr_.reset
    (
        new globalIndex(globalData.globalIndices())
    );
    const globalIndex& globalIndices = globalPointNumberingPtr_();

    // Create master to slave addressing. Empty for slave points.
    globalPointSlavesPtr_.reset
    (
        new labelListList(coupledPatch().nPoints())
    );
    labelListList& globalPointSlaves = globalPointSlavesPtr_();

    forAllConstIter(Map<label>, meshToProcPoint, iter)
    {
        label localPointI = iter.key();
        const labelList& pPoints = globalData.procPoints()[iter()];

        // Am I master?
        if
        (
            globalIndices.isLocal(pPoints[0])
         && globalIndices.toLocal(pPoints[0]) == localPointI
        )
        {
            labelList& slaves = globalPointSlaves[localPointI];
            slaves.setSize(pPoints.size()-1);
            for (label i = 1; i < pPoints.size(); i++)
            {
                slaves[i-1] = pPoints[i];
            }
        }
    }

    // Create schedule to get information from slaves onto master

    // Construct compact numbering and distribution map.
    // Changes globalPointSlaves to be indices into compact data

    List<Map<label> > compactMap(Pstream::nProcs());
    globalPointSlavesMapPtr_.reset
    (
        new mapDistribute
        (
            globalIndices,
            globalPointSlaves,
            compactMap
        )
    );

    if (debug)
    {
        Pout<< "globalMeshData::calcGlobalPointSlaves() :"
            << " coupled points:" << coupledPatch().nPoints()
            << " additional remote points:"
            <<  globalPointSlavesMapPtr_().constructSize()
              - coupledPatch().nPoints()
            << endl;
    }
}


void Foam::globalMeshData::calcGlobalEdgeSlaves() const
{
    if (debug)
    {
        Pout<< "globalMeshData::calcGlobalEdgeSlaves() :"
            << " calculating coupled master to slave edge addressing."
            << endl;
    }

    const labelListList& globalPointSlaves = this->globalPointSlaves();
    const mapDistribute& globalPointSlavesMap = this->globalPointSlavesMap();

    // - Send across connected edges (in global edge addressing)
    // - Check on receiving side whether edge has same slave edge
    //   on both endpoints.

    // Create global numbering for coupled edges
    globalEdgeNumberingPtr_.reset
    (
        new globalIndex(coupledPatch().nEdges())
    );
    const globalIndex& globalIndices = globalEdgeNumberingPtr_();

    // Coupled point to global coupled edges.
    labelListList globalPointEdges(globalPointSlavesMap.constructSize());

    // Create local version
    const labelListList& pointEdges = coupledPatch().pointEdges();
    forAll(pointEdges, pointI)
    {
        const labelList& pEdges = pointEdges[pointI];
        labelList& globalPEdges = globalPointEdges[pointI];
        globalPEdges.setSize(pEdges.size());
        forAll(pEdges, i)
        {
            globalPEdges[i] = globalIndices.toGlobal(pEdges[i]);
        }
    }

    // Pull slave data to master
    globalPointSlavesMap.distribute(globalPointEdges);

    // Now check on master if any of my edges are also on slave.
    // This assumes that if slaves have a coupled edge it is also on
    // the master (otherwise the mesh would be illegal anyway)

    labelHashSet pointEdgeSet;

    const edgeList& edges = coupledPatch().edges();

    // Create master to slave addressing. Empty for slave edges.
    globalEdgeSlavesPtr_.reset(new labelListList(edges.size()));
    labelListList& globalEdgeSlaves = globalEdgeSlavesPtr_();

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];
        const labelList& slaves0 = globalPointSlaves[e[0]];
        const labelList& slaves1 = globalPointSlaves[e[1]];

        // Check for edges that are in both slaves0 and slaves1.
        pointEdgeSet.clear();
        forAll(slaves0, i)
        {
            const labelList& connectedEdges = globalPointEdges[slaves0[i]];
            pointEdgeSet.insert(connectedEdges);
        }
        forAll(slaves1, i)
        {
            const labelList& connectedEdges = globalPointEdges[slaves1[i]];

            forAll(connectedEdges, j)
            {
                label globalEdgeI = connectedEdges[j];

                if (pointEdgeSet.found(globalEdgeI))
                {
                    // Found slave edge.
                    label sz = globalEdgeSlaves[edgeI].size();
                    globalEdgeSlaves[edgeI].setSize(sz+1);
                    globalEdgeSlaves[edgeI][sz] = globalEdgeI;
                }
            }
        }
    }


    // Construct map
    List<Map<label> > compactMap(Pstream::nProcs());
    globalEdgeSlavesMapPtr_.reset
    (
        new mapDistribute
        (
            globalIndices,
            globalEdgeSlaves,
            compactMap
        )
    );

    if (debug)
    {
        Pout<< "globalMeshData::calcGlobalEdgeSlaves() :"
            << " coupled edge:" << edges.size()
            << " additional remote edges:"
            << globalEdgeSlavesMapPtr_().constructSize() - edges.size()
            << endl;
    }
}


// Calculate uncoupled boundary faces (without calculating
// primitiveMesh::pointFaces())
void Foam::globalMeshData::calcPointBoundaryFaces
(
    labelListList& pointBoundaryFaces
) const
{
    const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();
    const Map<label>& meshPointMap = coupledPatch().meshPointMap();

    // 1. Count

    labelList nPointFaces(coupledPatch().nPoints(), 0);

    forAll(bMesh, patchI)
    {
        const polyPatch& pp = bMesh[patchI];

        if (!pp.coupled())
        {
            forAll(pp, i)
            {
                const face& f = pp[i];

                forAll(f, fp)
                {
                    Map<label>::const_iterator iter = meshPointMap.find(f[fp]);
                    if (iter != meshPointMap.end())
                    {
                        nPointFaces[iter()]++;
                    }
                }
            }
        }
    }


    // 2. Size

    pointBoundaryFaces.setSize(coupledPatch().nPoints());
    forAll(nPointFaces, pointI)
    {
        pointBoundaryFaces[pointI].setSize(nPointFaces[pointI]);
    }
    nPointFaces = 0;


    // 3. Fill

    forAll(bMesh, patchI)
    {
        const polyPatch& pp = bMesh[patchI];

        if (!pp.coupled())
        {
            forAll(pp, i)
            {
                const face& f = pp[i];
                forAll(f, fp)
                {
                    Map<label>::const_iterator iter = meshPointMap.find(f[fp]);
                    if (iter != meshPointMap.end())
                    {
                        label bFaceI = pp.start() + i - mesh_.nInternalFaces();
                        pointBoundaryFaces[iter()][nPointFaces[iter()]++] =
                            bFaceI;
                    }
                }
            }
        }
    }
}

void Foam::globalMeshData::calcGlobalPointBoundaryFaces() const
{
    if (debug)
    {
        Pout<< "globalMeshData::calcGlobalPointBoundaryFaces() :"
            << " calculating coupled point to boundary face addressing."
            << endl;
    }

    // Construct local point to (uncoupled)boundaryfaces.
    labelListList pointBoundaryFaces;
    calcPointBoundaryFaces(pointBoundaryFaces);


    // Global indices for boundary faces
    globalBoundaryFaceNumberingPtr_.reset
    (
        new globalIndex(mesh_.nFaces()-mesh_.nInternalFaces())
    );
    globalIndex& globalIndices = globalBoundaryFaceNumberingPtr_();


    // Convert local boundary faces to global numbering
    globalPointBoundaryFacesPtr_.reset
    (
        new labelListList(globalPointSlavesMap().constructSize())
    );
    labelListList& globalPointBoundaryFaces = globalPointBoundaryFacesPtr_();

    forAll(pointBoundaryFaces, pointI)
    {
        const labelList& bFaces = pointBoundaryFaces[pointI];
        labelList& globalFaces = globalPointBoundaryFaces[pointI];
        globalFaces.setSize(bFaces.size());
        forAll(bFaces, i)
        {
            globalFaces[i] = globalIndices.toGlobal(bFaces[i]);
        }
    }


    // Pull slave pointBoundaryFaces to master
    globalPointSlavesMap().distribute(globalPointBoundaryFaces);


    // Merge slave labels into master globalPointBoundaryFaces
    const labelListList& pointSlaves = globalPointSlaves();

    forAll(pointSlaves, pointI)
    {
        const labelList& slaves = pointSlaves[pointI];

        if (slaves.size() > 0)
        {
            labelList& myBFaces = globalPointBoundaryFaces[pointI];

            forAll(slaves, i)
            {
                const labelList& slaveBFaces =
                    globalPointBoundaryFaces[slaves[i]];

                // Add all slaveBFaces. Note that need to check for
                // uniqueness only in case of cyclics.

                label sz = myBFaces.size();
                myBFaces.setSize(sz+slaveBFaces.size());
                forAll(slaveBFaces, j)
                {
                    label slave = slaveBFaces[j];
                    if (findIndex(SubList<label>(myBFaces, sz), slave) == -1)
                    {
                        myBFaces[sz++] = slave;
                    }
                }
                myBFaces.setSize(sz);
            }
        }
    }


    // Copy merged boundaryFaces back from master into slave slot
    forAll(pointSlaves, pointI)
    {
        const labelList& bFaces = globalPointBoundaryFaces[pointI];
        const labelList& slaves = pointSlaves[pointI];

        forAll(slaves, i)
        {
            globalPointBoundaryFaces[slaves[i]] = bFaces;
        }
    }


    // Sync back to slaves.
    globalPointSlavesMap().reverseDistribute
    (
        coupledPatch().nPoints(),
        globalPointBoundaryFaces
    );


    // Construct a map to get the face data directly
    List<Map<label> > compactMap(Pstream::nProcs());
    globalPointBoundaryFacesMapPtr_.reset
    (
        new mapDistribute
        (
            globalIndices,
            globalPointBoundaryFaces,
            compactMap
        )
    );

    if (debug)
    {
        Pout<< "globalMeshData::calcGlobalPointBoundaryFaces() :"
            << " coupled points:" << coupledPatch().nPoints()
            << " local boundary faces:" <<  globalIndices.localSize()
            << " additional remote faces:"
            <<  globalPointBoundaryFacesMapPtr_().constructSize()
              - globalIndices.localSize()
            << endl;
    }
}


void Foam::globalMeshData::calcGlobalPointBoundaryCells() const
{
    if (debug)
    {
        Pout<< "globalMeshData::calcGlobalPointBoundaryCells() :"
            << " calculating coupled point to boundary cell addressing."
            << endl;
    }

    // Create map of boundary cells and point-cell addressing
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    label bCellI = 0;
    Map<label> meshCellMap(4*coupledPatch().nPoints());
    DynamicList<label> cellMap(meshCellMap.size());

    // Create addressing for point to boundary cells (local)
    labelListList pointBoundaryCells(coupledPatch().nPoints());

    forAll(coupledPatch().meshPoints(), pointI)
    {
        label meshPointI = coupledPatch().meshPoints()[pointI];
        const labelList& pCells = mesh_.pointCells(meshPointI);

        labelList& bCells = pointBoundaryCells[pointI];
        bCells.setSize(pCells.size());

        forAll(pCells, i)
        {
            label cellI = pCells[i];
            Map<label>::iterator fnd = meshCellMap.find(cellI);

            if (fnd != meshCellMap.end())
            {
                bCells[i] = fnd();
            }
            else
            {
                meshCellMap.insert(cellI, bCellI);
                cellMap.append(cellI);
                bCells[i] = bCellI;
                bCellI++;
            }
        }
    }


    boundaryCellsPtr_.reset(new labelList());
    labelList& boundaryCells = boundaryCellsPtr_();
    boundaryCells.transfer(cellMap.shrink());


    // Convert point-cells to global point numbers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    globalBoundaryCellNumberingPtr_.reset
    (
        new globalIndex(boundaryCells.size())
    );
    globalIndex& globalIndices = globalBoundaryCellNumberingPtr_();


    globalPointBoundaryCellsPtr_.reset
    (
        new labelListList(globalPointSlavesMap().constructSize())
    );
    labelListList& globalPointBoundaryCells = globalPointBoundaryCellsPtr_();

    forAll(pointBoundaryCells, pointI)
    {
        const labelList& pCells = pointBoundaryCells[pointI];
        labelList& globalCells = globalPointBoundaryCells[pointI];
        globalCells.setSize(pCells.size());
        forAll(pCells, i)
        {
            globalCells[i] = globalIndices.toGlobal(pCells[i]);
        }
    }


    // Pull slave pointBoundaryCells to master
    globalPointSlavesMap().distribute(globalPointBoundaryCells);


    // Merge slave labels into master globalPointBoundaryCells
    const labelListList& pointSlaves = globalPointSlaves();

    forAll(pointSlaves, pointI)
    {
        const labelList& slaves = pointSlaves[pointI];

        if (slaves.size() > 0)
        {
            labelList& myBCells = globalPointBoundaryCells[pointI];

            forAll(slaves, i)
            {
                const labelList& slaveBCells =
                    globalPointBoundaryCells[slaves[i]];

                // Add all slaveBCells. Note that need to check for
                // uniqueness only in case of cyclics.

                label sz = myBCells.size();
                myBCells.setSize(sz+slaveBCells.size());
                forAll(slaveBCells, j)
                {
                    label slave = slaveBCells[j];
                    if (findIndex(SubList<label>(myBCells, sz), slave) == -1)
                    {
                        myBCells[sz++] = slave;
                    }
                }
                myBCells.setSize(sz);
            }
        }
    }


    // Copy merged boundaryCells back from master into slave slot
    forAll(pointSlaves, pointI)
    {
        const labelList& bCells = globalPointBoundaryCells[pointI];
        const labelList& slaves = pointSlaves[pointI];

        forAll(slaves, i)
        {
            globalPointBoundaryCells[slaves[i]] = bCells;
        }
    }


    // Sync back to slaves.
    globalPointSlavesMap().reverseDistribute
    (
        coupledPatch().nPoints(),
        globalPointBoundaryCells
    );


    // Construct a map to get the cell data directly
    List<Map<label> > compactMap(Pstream::nProcs());
    globalPointBoundaryCellsMapPtr_.reset
    (
        new mapDistribute
        (
            globalIndices,
            globalPointBoundaryCells,
            compactMap
        )
    );

    if (debug)
    {
        Pout<< "globalMeshData::calcGlobalPointBoundaryCells() :"
            << " coupled points:" << coupledPatch().nPoints()
            << " local boundary cells:" <<  globalIndices.localSize()
            << " additional remote cells:"
            <<  globalPointBoundaryCellsMapPtr_().constructSize()
              - globalIndices.localSize()
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from polyMesh
Foam::globalMeshData::globalMeshData(const polyMesh& mesh)
:
    processorTopology(mesh.boundaryMesh()),
    mesh_(mesh),
    bb_(vector::zero, vector::zero),
    nTotalPoints_(-1),
    nTotalFaces_(-1),
    nTotalCells_(-1),
    processorPatches_(0),
    processorPatchIndices_(0),
    processorPatchNeighbours_(0),
    nGlobalPoints_(-1),
    sharedPointLabels_(0),
    sharedPointAddr_(0),
    sharedPointGlobalLabelsPtr_(NULL),
    nGlobalEdges_(-1),
    sharedEdgeLabelsPtr_(NULL),
    sharedEdgeAddrPtr_(NULL)
{
    updateMesh();
}


// Read constructor given IOobject and a polyMesh reference
Foam::globalMeshData::globalMeshData(const IOobject& io, const polyMesh& mesh)
:
    processorTopology(mesh.boundaryMesh()),
    mesh_(mesh),
    bb_(mesh.points()),
    nTotalPoints_(-1),
    nTotalFaces_(-1),
    nTotalCells_(-1),
    processorPatches_(0),
    processorPatchIndices_(0),
    processorPatchNeighbours_(0),
    nGlobalPoints_(-1),
    sharedPointLabels_(0),
    sharedPointAddr_(0),
    sharedPointGlobalLabelsPtr_(NULL),
    nGlobalEdges_(-1),
    sharedEdgeLabelsPtr_(NULL),
    sharedEdgeAddrPtr_(NULL)
{
    initProcAddr();

    IOdictionary dict(io);

    dict.lookup("nTotalPoints") >> nTotalPoints_;
    dict.lookup("nTotalFaces") >> nTotalFaces_;
    dict.lookup("nTotalCells") >> nTotalCells_;
    dict.lookup("nGlobalPoints") >> nGlobalPoints_;
    dict.lookup("sharedPointLabels") >> sharedPointLabels_;
    dict.lookup("sharedPointAddr") >> sharedPointAddr_;
    labelList sharedPointGlobalLabels(dict.lookup("sharedPointGlobalLabels"));

    sharedPointGlobalLabelsPtr_ = new labelList(sharedPointGlobalLabels);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::globalMeshData::~globalMeshData()
{
    clearOut();
}


void Foam::globalMeshData::clearOut()
{
    deleteDemandDrivenData(sharedPointGlobalLabelsPtr_);
    // Edge
    nGlobalPoints_ = -1;
    deleteDemandDrivenData(sharedEdgeLabelsPtr_);
    deleteDemandDrivenData(sharedEdgeAddrPtr_);

    coupledPatchPtr_.clear();
    // Point
    globalPointNumberingPtr_.clear();
    globalPointSlavesPtr_.clear();
    globalPointSlavesMapPtr_.clear();
    // Edge
    globalEdgeNumberingPtr_.clear();
    globalEdgeSlavesPtr_.clear();
    globalEdgeSlavesMapPtr_.clear();
    // Face
    globalBoundaryFaceNumberingPtr_.clear();
    globalPointBoundaryFacesPtr_.clear();
    globalPointBoundaryFacesMapPtr_.clear();
    // Cell
    boundaryCellsPtr_.clear();
    globalBoundaryCellNumberingPtr_.clear();
    globalPointBoundaryCellsPtr_.clear();
    globalPointBoundaryCellsMapPtr_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return shared point global labels.
const Foam::labelList& Foam::globalMeshData::sharedPointGlobalLabels() const
{
    if (!sharedPointGlobalLabelsPtr_)
    {
        sharedPointGlobalLabelsPtr_ = new labelList(sharedPointLabels_.size());
        labelList& sharedPointGlobalLabels = *sharedPointGlobalLabelsPtr_;

        IOobject addrHeader
        (
            "pointProcAddressing",
            mesh_.facesInstance()/mesh_.meshSubDir,
            mesh_,
            IOobject::MUST_READ
        );

        if (addrHeader.headerOk())
        {
            // There is a pointProcAddressing file so use it to get labels
            // on the original mesh
            Pout<< "globalMeshData::sharedPointGlobalLabels : "
                << "Reading pointProcAddressing" << endl;

            labelIOList pointProcAddressing(addrHeader);

            forAll(sharedPointLabels_, i)
            {
                // Get my mesh point
                label pointI = sharedPointLabels_[i];

                // Map to mesh point of original mesh
                sharedPointGlobalLabels[i] = pointProcAddressing[pointI];
            }
        }
        else
        {
            Pout<< "globalMeshData::sharedPointGlobalLabels :"
                << " Setting pointProcAddressing to -1" << endl;

            sharedPointGlobalLabels = -1;
        }
    }
    return *sharedPointGlobalLabelsPtr_;
}


// Collect coordinates of shared points. (does parallel communication!)
Foam::pointField Foam::globalMeshData::sharedPoints() const
{
    // Get all processors to send their shared points to master.
    // (not very efficient)

    pointField sharedPoints(nGlobalPoints_);

    if (Pstream::master())
    {
        // Master:
        // insert my own data first
        forAll(sharedPointLabels_, i)
        {
            label sharedPointI = sharedPointAddr_[i];

            sharedPoints[sharedPointI] = mesh_.points()[sharedPointLabels_[i]];
        }

        // Receive data from slaves and insert
        for
        (
            int slave=Pstream::firstSlave();
            slave<=Pstream::lastSlave();
            slave++
        )
        {
            IPstream fromSlave(Pstream::blocking, slave);

            labelList nbrSharedPointAddr;
            pointField nbrSharedPoints;
            fromSlave >> nbrSharedPointAddr >> nbrSharedPoints;

            forAll(nbrSharedPointAddr, i)
            {
                label sharedPointI = nbrSharedPointAddr[i];

                sharedPoints[sharedPointI] = nbrSharedPoints[i];
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
            OPstream toSlave
            (
                Pstream::blocking,
                slave,
                sharedPoints.size()*sizeof(vector::zero)
            );
            toSlave << sharedPoints;
        }
    }
    else
    {
        // Slave:
        // send points
        {
            OPstream toMaster(Pstream::blocking, Pstream::masterNo());

            toMaster
                << sharedPointAddr_
                << UIndirectList<point>(mesh_.points(), sharedPointLabels_)();
        }

        // Receive sharedPoints
        {
            IPstream fromMaster(Pstream::blocking, Pstream::masterNo());
            fromMaster >> sharedPoints;
        }
    }

    return sharedPoints;
}


// Collect coordinates of shared points. (does parallel communication!)
Foam::pointField Foam::globalMeshData::geometricSharedPoints() const
{
    // Get coords of my shared points
    pointField sharedPoints(sharedPointLabels_.size());

    forAll(sharedPointLabels_, i)
    {
        label meshPointI = sharedPointLabels_[i];

        sharedPoints[i] = mesh_.points()[meshPointI];
    }

    // Append from all processors
    combineReduce(sharedPoints, plusEqOp<pointField>());

    // Merge tolerance
    scalar tolDim = matchTol_ * bb_.mag();

    // And see how many are unique
    labelList pMap;
    pointField mergedPoints;

    mergePoints
    (
        sharedPoints,   // coordinates to merge
        tolDim,         // tolerance
        false,          // verbosity
        pMap,
        mergedPoints
    );

    return mergedPoints;
}


Foam::label Foam::globalMeshData::nGlobalEdges() const
{
    if (nGlobalEdges_ == -1)
    {
        calcSharedEdges();
    }
    return nGlobalEdges_;
}


const Foam::labelList& Foam::globalMeshData::sharedEdgeLabels() const
{
    if (!sharedEdgeLabelsPtr_)
    {
        calcSharedEdges();
    }
    return *sharedEdgeLabelsPtr_;
}


const Foam::labelList& Foam::globalMeshData::sharedEdgeAddr() const
{
    if (!sharedEdgeAddrPtr_)
    {
        calcSharedEdges();
    }
    return *sharedEdgeAddrPtr_;
}


const Foam::indirectPrimitivePatch& Foam::globalMeshData::coupledPatch() const
{
    if (!coupledPatchPtr_.valid())
    {
        const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();

        label nCoupled = 0;

        forAll(bMesh, patchI)
        {
            const polyPatch& pp = bMesh[patchI];

            if (pp.coupled())
            {
                nCoupled += pp.size();
            }
        }
        labelList coupledFaces(nCoupled);
        nCoupled = 0;

        forAll(bMesh, patchI)
        {
            const polyPatch& pp = bMesh[patchI];

            if (pp.coupled())
            {
                label faceI = pp.start();

                forAll(pp, i)
                {
                    coupledFaces[nCoupled++] = faceI++;
                }
            }
        }

        coupledPatchPtr_.reset
        (
            new indirectPrimitivePatch
            (
                IndirectList<face>
                (
                    mesh_.faces(),
                    coupledFaces
                ),
                mesh_.points()
            )
        );

        if (debug)
        {
            Pout<< "globalMeshData::coupledPatch() :"
                << " constructed  coupled faces patch:"
                << "  faces:" << coupledPatchPtr_().size()
                << "  points:" << coupledPatchPtr_().nPoints()
                << endl;
        }
    }
    return coupledPatchPtr_();
}


const Foam::globalIndex& Foam::globalMeshData::globalPointNumbering() const
{
    if (!globalPointNumberingPtr_.valid())
    {
        calcGlobalPointSlaves();
    }
    return globalPointNumberingPtr_();
}


const Foam::labelListList& Foam::globalMeshData::globalPointSlaves() const
{
    if (!globalPointSlavesPtr_.valid())
    {
        calcGlobalPointSlaves();
    }
    return globalPointSlavesPtr_();
}


const Foam::mapDistribute& Foam::globalMeshData::globalPointSlavesMap() const
{
    if (!globalPointSlavesMapPtr_.valid())
    {
        calcGlobalPointSlaves();
    }
    return globalPointSlavesMapPtr_();
}


const Foam::globalIndex& Foam::globalMeshData::globalEdgeNumbering() const
{
    if (!globalEdgeNumberingPtr_.valid())
    {
        calcGlobalEdgeSlaves();
    }
    return globalEdgeNumberingPtr_();
}


const Foam::labelListList& Foam::globalMeshData::globalEdgeSlaves() const
{
    if (!globalEdgeSlavesPtr_.valid())
    {
        calcGlobalEdgeSlaves();
    }
    return globalEdgeSlavesPtr_();
}


const Foam::mapDistribute& Foam::globalMeshData::globalEdgeSlavesMap() const
{
    if (!globalEdgeSlavesMapPtr_.valid())
    {
        calcGlobalEdgeSlaves();
    }
    return globalEdgeSlavesMapPtr_();
}


const Foam::globalIndex& Foam::globalMeshData::globalBoundaryFaceNumbering()
const
{
    if (!globalBoundaryFaceNumberingPtr_.valid())
    {
        calcGlobalPointBoundaryFaces();
    }
    return globalBoundaryFaceNumberingPtr_();
}


const Foam::labelListList& Foam::globalMeshData::globalPointBoundaryFaces()
const
{
    if (!globalPointBoundaryFacesPtr_.valid())
    {
        calcGlobalPointBoundaryFaces();
    }
    return globalPointBoundaryFacesPtr_();
}


const Foam::mapDistribute& Foam::globalMeshData::globalPointBoundaryFacesMap()
const
{
    if (!globalPointBoundaryFacesMapPtr_.valid())
    {
        calcGlobalPointBoundaryFaces();
    }
    return globalPointBoundaryFacesMapPtr_();
}


const Foam::labelList& Foam::globalMeshData::boundaryCells() const
{
    if (!boundaryCellsPtr_.valid())
    {
        calcGlobalPointBoundaryCells();
    }
    return boundaryCellsPtr_();
}


const Foam::globalIndex& Foam::globalMeshData::globalBoundaryCellNumbering()
const
{
    if (!globalBoundaryCellNumberingPtr_.valid())
    {
        calcGlobalPointBoundaryCells();
    }
    return globalBoundaryCellNumberingPtr_();
}


const Foam::labelListList& Foam::globalMeshData::globalPointBoundaryCells()
const
{
    if (!globalPointBoundaryCellsPtr_.valid())
    {
        calcGlobalPointBoundaryCells();
    }
    return globalPointBoundaryCellsPtr_();
}


const Foam::mapDistribute& Foam::globalMeshData::globalPointBoundaryCellsMap()
const
{
    if (!globalPointBoundaryCellsMapPtr_.valid())
    {
        calcGlobalPointBoundaryCells();
    }
    return globalPointBoundaryCellsMapPtr_();
}


void Foam::globalMeshData::movePoints(const pointField& newPoints)
{
    // Topology does not change and we don't store any geometry so nothing
    // needs to be done.
}


// Update all data after morph
void Foam::globalMeshData::updateMesh()
{
    // Clear out old data
    clearOut();

    // Do processor patch addressing
    initProcAddr();

    // Note: boundBox does reduce
    bb_ = boundBox(mesh_.points());

    scalar tolDim = matchTol_ * bb_.mag();

    if (debug)
    {
        Pout<< "globalMeshData : bb_:" << bb_
            << " merge dist:" << tolDim << endl;
    }


    // Option 1. Topological
    {
        // Calculate all shared points. This does all the hard work.
        globalPoints parallelPoints(mesh_, false);

        // Copy data out.
        nGlobalPoints_ = parallelPoints.nGlobalPoints();
        sharedPointLabels_ = parallelPoints.sharedPointLabels();
        sharedPointAddr_ = parallelPoints.sharedPointAddr();
    }
    //// Option 2. Geometric
    //{
    //    // Calculate all shared points. This does all the hard work.
    //    geomGlobalPoints parallelPoints(mesh_, tolDim);
    //
    //    // Copy data out.
    //    nGlobalPoints_ = parallelPoints.nGlobalPoints();
    //    sharedPointLabels_ = parallelPoints.sharedPointLabels();
    //    sharedPointAddr_ = parallelPoints.sharedPointAddr();
    //
    //    nGlobalEdges_ = parallelPoints.nGlobalEdges();
    //    sharedEdgeLabels_ = parallelPoints.sharedEdgeLabels();
    //    sharedEdgeAddr_ = parallelPoints.sharedEdgeAddr();
    //}

    // Total number of faces. Start off from all faces. Remove coincident
    // processor faces (on highest numbered processor) before summing.
    nTotalFaces_ = mesh_.nFaces();

    // Do not count processor-patch faces that are coincident.
    forAll(processorPatches_, i)
    {
        label patchI = processorPatches_[i];

        if (isType<processorPolyPatch>(mesh_.boundaryMesh()[patchI]))
        {
            // Normal, unseparated processor patch. Remove duplicates.
            nTotalFaces_ -= mesh_.boundaryMesh()[patchI].size();
        }
    }
    reduce(nTotalFaces_, sumOp<label>());

    if (debug)
    {
        Pout<< "globalMeshData : nTotalFaces_:" << nTotalFaces_ << endl;
    }


    nTotalCells_ = mesh_.nCells();
    reduce(nTotalCells_, sumOp<label>());

    if (debug)
    {
        Pout<< "globalMeshData : nTotalCells_:" << nTotalCells_ << endl;
    }

    nTotalPoints_ = mesh_.nPoints();

    // Correct points for duplicate ones. We have
    // - points shared between 2 processor patches only. Count only on
    //   lower numbered processor. Make sure to count only once since points
    //   can be on multiple patches on the same processor.
    // - globally shared points.

    if (Pstream::parRun())
    {
        const label UNSET = 0;      // not set
        const label SHARED = 1;     // globally shared
        const label VISITED = 2;    // corrected for

        // Mark globally shared points
        PackedList<2> pointStatus(mesh_.nPoints(), UNSET);

        forAll(sharedPointLabels_, i)
        {
            label meshPointI = sharedPointLabels_[i];

            pointStatus.set(meshPointI, SHARED);
        }

        // Send patch local points
        forAll(processorPatches_, i)
        {
            label patchI = processorPatches_[i];

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(mesh_.boundaryMesh()[patchI]);

            OPstream toNeighbour(Pstream::blocking, procPatch.neighbProcNo());

            toNeighbour << procPatch.localPoints();
        }

        // Receive patch local points and uncount if coincident (and not shared)
        forAll(processorPatches_, i)
        {
            label patchI = processorPatches_[i];

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(mesh_.boundaryMesh()[patchI]);

            IPstream fromNeighbour(Pstream::blocking, procPatch.neighbProcNo());

            pointField nbrPoints(fromNeighbour);

            if (Pstream::myProcNo() > procPatch.neighbProcNo())
            {
                labelList pMap;
                matchPoints
                (
                    procPatch.localPoints(),
                    nbrPoints,
                    scalarField(procPatch.nPoints(), tolDim),   // tolerance
                    false,      // verbosity
                    pMap        // map from my points to nbrPoints
                );

                forAll(pMap, patchPointI)
                {
                    label meshPointI = procPatch.meshPoints()[patchPointI];

                    label stat = pointStatus.get(meshPointI);

                    if (stat == UNSET)
                    {
                        // Mark point as visited so if point is on multiple proc
                        // patches it only gets uncounted once.
                        pointStatus.set(meshPointI, VISITED);

                        if (pMap[patchPointI] != -1)
                        {
                            // Points share same coordinate so uncount.
                            nTotalPoints_--;
                        }
                    }
                }
            }
        }
        // Sum all points
        reduce(nTotalPoints_, sumOp<label>());
    }

    // nTotalPoints has not been corrected yet for shared points. For these
    // just collect all their coordinates and count unique ones.

    label mySharedPoints = sharedPointLabels_.size();
    reduce(mySharedPoints, sumOp<label>());

    // Collect and merge shared points (does parallel communication)
    pointField geomSharedPoints(geometricSharedPoints());
    label nGeomSharedPoints = geomSharedPoints.size();

    // Shared points merged down to mergedPoints size.
    nTotalPoints_ -= mySharedPoints - nGeomSharedPoints;

    if (debug)
    {
        Pout<< "globalMeshData : nTotalPoints_:" << nTotalPoints_ << endl;
    }

    //
    // Now we have all info we wanted.
    // Do some checking (if debug is set)
    //

    if (debug)
    {
        if (Pstream::master())
        {
            // We have the geometricSharedPoints already so write them.
            // Ideally would like to write the networks of connected points as
            // well but this is harder. (Todo)
            Pout<< "globalMeshData : writing geometrically separated shared"
                << " points to geomSharedPoints.obj" << endl;

            OFstream str("geomSharedPoints.obj");

            forAll(geomSharedPoints, i)
            {
                const point& pt = geomSharedPoints[i];

                str << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z()
                    << nl;
            }
        }
    }
}


// Write data
bool Foam::globalMeshData::write() const
{
    IOdictionary dict
    (
        IOobject
        (
            "parallelData",
            mesh_.facesInstance(),
            mesh_.meshSubDir,
            mesh_
        )
    );

    dict.add("nTotalPoints", nTotalPoints());
    dict.add("nTotalFaces", nTotalFaces());
    dict.add("nTotalCells", nTotalCells());

    dict.add("nGlobalPoints", nGlobalPoints());
    dict.add("sharedPointLabels", sharedPointLabels());
    dict.add("sharedPointAddr", sharedPointAddr());
    dict.add("sharedPointGlobalLabels", sharedPointGlobalLabels());

    return dict.writeObject
    (
        IOstream::ASCII,
        IOstream::currentVersion,
        IOstream::UNCOMPRESSED
    );
}


// * * * * * * * * * * * * * * * Ostream Operators * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const globalMeshData& p)
{
    os  << "nTotalPoints " << p.nTotalPoints() << token::END_STATEMENT << nl
        << "nTotalFaces " << p.nTotalFaces() << token::END_STATEMENT << nl
        << "nTotalCells " << p.nTotalCells() << token::END_STATEMENT << nl
        << "nGlobalPoints " << p.nGlobalPoints() << token::END_STATEMENT << nl
        << "sharedPointLabels " << p.sharedPointLabels()
        << token::END_STATEMENT << nl
        << "sharedPointAddr " << p.sharedPointAddr()
        << token::END_STATEMENT << endl;

    return os;
}


// ************************************************************************* //
