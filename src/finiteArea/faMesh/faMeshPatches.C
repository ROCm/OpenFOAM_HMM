/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 Wikki Ltd
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "faMesh.H"
#include "IndirectList.H"
#include "faPatchData.H"
#include "processorPolyPatch.H"
#include "processorFaPatch.H"
#include "globalMeshData.H"
#include "indirectPrimitivePatch.H"
#include "edgeHashes.H"
#include "LabelledItem.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Manage patch pairs with a 'labelled' edge.
// The edge first/second correspond to the owner/neighbour patches.
// The index is a face index on the neighbour patch (FUTURE).

// Local typedefs

typedef LabelledItem<edge> patchPairInfo;
typedef List<patchPairInfo> patchPairInfoList;
typedef UIndirectList<patchPairInfo> patchPairInfoUIndList;


// Handling of dangling coupled edges.
// Tag values to "push" with special -(patchId+2)
struct combineDanglingEdge
{
    const label upperLimit;

    // Set dangling patchId from real patchId
    static void setDangling(patchPairInfo& pairing, const label patchId)
    {
        pairing.first() = pairing.second() = -(patchId + 2);
        pairing.setIndex(-1);  // Invalidate
    }

    // Convert dangling patchId to real patchId
    static void correct(patchPairInfo& pairing)
    {
        if (pairing.first() < -1)
        {
            pairing.first() = -(pairing.first() + 2);
        }
        if (pairing.second() < -1)
        {
            pairing.second() = -(pairing.second() + 2);
        }
    }

    //- Construct with upper limit (the number of non-processor patches)
    explicit combineDanglingEdge(const label nNonProcessor)
    :
        upperLimit(nNonProcessor)
    {}


    // Combine operation: overwrite unused or processor patches with
    // 'dangling' patch information only
    void operator()(patchPairInfo& x, const patchPairInfo& y) const
    {
        if (y.first() < -1 && edge::compare(x, y) == 0)
        {
            if (x.first() == -1 || x.first() >= upperLimit)
            {
                x.first() = y.first();
            }
            if (x.second() == -1 || x.second() >= upperLimit)
            {
                x.second() = y.first();
                x.index() = y.index();
            }
        }
    }
};


// Populate patch pairings according to the boundary edges
void findEdgePatchPairing
(
    const polyBoundaryMesh& pbm,
    const EdgeMap<label>& edgeToIndex,
    patchPairInfoList& patchPairs,
    label nMissing = -1
)
{
    // Count how many slots (both sides) to be filled
    if (nMissing < 0)
    {
        nMissing = 0;
        for (const patchPairInfo& pairing : patchPairs)
        {
            if (pairing.first() == -1) ++nMissing;
            if (pairing.second() == -1) ++nMissing;
        }
    }

    forAll(pbm, patchID)
    {
        if (!nMissing) break;  // Everything filled

        const polyPatch& pp = pbm[patchID];

        const bool isProcPatch = isA<processorPolyPatch>(pp);

        // Examine neighbour boundary edges
        for (label edgei = pp.nInternalEdges(); edgei < pp.nEdges(); ++edgei)
        {
            // Lookup global edge of this neighbour,
            // find matching owner boundary edge

            const label edgeIndex =
                edgeToIndex.lookup(pp.meshEdge(edgei), -1);

            if (edgeIndex != -1)
            {
                // Add patchId.
                // - hash-like so will only insert once
                // - also add in the attached patch face (there is only one),
                //   saved in mesh face numbering for processor patches

                patchPairInfo& pairing = patchPairs[edgeIndex];

                if (pairing.insert(patchID))
                {
                    if (isProcPatch)
                    {
                        // Save the mesh face

                        const label patchFacei = pp.edgeFaces()[edgei][0];

                        pairing.index() = (pp.start() + patchFacei);
                    }

                    --nMissing;

                    if (!nMissing) break;  // Early exit
                }
            }
        }
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faMesh::reorderProcEdges
(
    faPatchData& patchDef,
    const List<LabelledItem<edge>>& bndEdgePatchPairs
) const
{
    if (!patchDef.coupled() || patchDef.edgeLabels_.empty())
    {
        return;
    }

    const polyBoundaryMesh& pbm = mesh().boundaryMesh();

    const label procPatchID = patchDef.neighPolyPatchId_;

    const auto* procPatch = isA<processorPolyPatch>(pbm[procPatchID]);

    if (!procPatch)
    {
        FatalErrorInFunction
            << "Internal addressing error. Patch " << procPatchID
            << " is not a processor patch" << nl
            << abort(FatalError);
    }

    // Reorder processor edges using order of neighbour processorPolyPatch

    const label nProcEdges = patchDef.edgeLabels_.size();

    labelList procFaces(nProcEdges, -1);

    forAll(procFaces, edgei)
    {
        const label bndEdgei =
            (patchDef.edgeLabels_[edgei] - patch().nInternalEdges());

        procFaces[edgei] = bndEdgePatchPairs[bndEdgei].index();
    }

    // Ascending proc-face numbering
    const labelList sortIndices(Foam::sortedOrder(procFaces));

    if (procFaces.size() && procFaces[sortIndices[0]] < 0)
    {
        FatalErrorInFunction
            << "Internal addressing error. Patch " << procPatchID
            << " with negative face" << nl
            << abort(FatalError);
    }


    const labelList& oldEdgeLabels = patchDef.edgeLabels_;
    labelList newEdgeLabels(oldEdgeLabels.size());

    // Most of the time, an individual proc-face will only be singly
    // attached to the finite-area patch. In rarer case, there could
    // multiple connections. For these cases, need to walk the face
    // edges - the direction depends on owner vs neighbour side.

    EdgeMap<label> multihit;

    for (label edgei = 0; edgei < nProcEdges; /*nil*/)
    {
        const label meshFacei = procFaces[sortIndices[edgei]];

        // Find all identical faces
        label endEdgei = edgei + 1;  // one beyond
        while
        (
            (endEdgei < nProcEdges)
         && (meshFacei == procFaces[sortIndices[endEdgei]])
        )
        {
            ++endEdgei;
        }

        if (edgei + 1 == endEdgei)
        {
            // Simplest case - a single connection

            newEdgeLabels[edgei] = oldEdgeLabels[sortIndices[edgei]];
        }
        else
        {
            multihit.clear();

            // Map from global edge to patch local edgeId
            for (label i = edgei; i < endEdgei; ++i)
            {
                const label patchEdgei = oldEdgeLabels[sortIndices[i]];

                // The edge in mesh numbering
                multihit.insert
                (
                    patch().meshEdge(patchEdgei),
                    patchEdgei
                );
            }

            if (multihit.size() != (endEdgei - edgei))
            {
                FatalErrorInFunction
                    << "Could only hash " << multihit.size()
                    << " edges from " << (endEdgei - edgei)
                    << " ... indicates a non-manifold connection" << nl
                    << multihit << nl
                    << abort(FatalError);
            }

            const face& f = mesh().faces()[meshFacei];

            forAll(f, fedgei)  // Note size() == nEdges()
            {
                edge e =
                (
                    patchDef.owner()
                  ? f.edge(fedgei)      // Forward walk
                  : f.rcEdge(fedgei)    // Reverse walk
                );

                auto iter = multihit.find(e);
                if (iter.found())
                {
                    newEdgeLabels[edgei++] = iter.val();
                    multihit.erase(iter);
                    if (multihit.empty())
                    {
                        break;
                    }
                }
            }

            if (edgei != endEdgei)
            {
                FatalErrorInFunction
                    << "Missed " << (edgei < endEdgei)
                    << " edges for face: " << meshFacei
                    << " ... indicates serious geometry issue" << nl
                    << multihit << nl
                    << abort(FatalError);
            }
            if (!multihit.empty())
            {
                FatalErrorInFunction
                    << "Missed edges for face: " << meshFacei
                    << " ... indicates serious geometry issue" << nl
                    << multihit << nl
                    << abort(FatalError);
            }
        }
        edgei = endEdgei;
    }

    patchDef.edgeLabels_.transfer(newEdgeLabels);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faMesh::addFaPatches
(
    PtrList<faPatch>& plist,
    const bool validBoundary
)
{
    if (!boundary().empty())
    {
        FatalErrorInFunction
            << "boundary already exists"
            << abort(FatalError);
    }

    globalMeshDataPtr_.reset(nullptr);

    boundary_.transfer(plist);

    setPrimitiveMeshData();

    if (validBoundary)
    {
        boundary_.checkDefinition();
    }
}


void Foam::faMesh::addFaPatches
(
    const List<faPatch*>& p,
    const bool validBoundary
)
{
    // Acquire ownership of the pointers
    PtrList<faPatch> plist(const_cast<List<faPatch*>&>(p));

    addFaPatches(plist, validBoundary);
}


Foam::PtrList<Foam::faPatch> Foam::faMesh::createOnePatch
(
    const word& patchName,
    const word& patchType
) const
{
    dictionary onePatchDict;
    if (!patchName.empty())
    {
        onePatchDict.add("name", patchName);
    }
    if (!patchType.empty())
    {
        onePatchDict.add("type", patchType);
    }

    return createPatchList
    (
        dictionary::null,
        word::null,     // Name for empty patch placeholder
        &onePatchDict   // Definitions for defaultPatch
    );
}


Foam::List<Foam::LabelledItem<Foam::edge>>
Foam::faMesh::getBoundaryEdgePatchPairs() const
{
    const polyBoundaryMesh& pbm = mesh().boundaryMesh();

    const label nInternalEdges = patch().nInternalEdges();
    const label nBoundaryEdges = patch().nBoundaryEdges();

    // Map edges (mesh numbering) back to a boundary index
    EdgeMap<label> edgeToBoundaryIndex(2*nBoundaryEdges);

    // Use labelled 'edge' for accounting of patch pairs
    patchPairInfoList bndEdgePatchPairs(nBoundaryEdges);


    // Pass 1:
    // - setup lookup (edge -> bnd index)
    // - add owner patch for each boundary edge
    {
        const SubList<labelList> bndEdgeToPatchFace
        (
            patch().edgeFaces(),
            patch().nBoundaryEdges(),
            patch().nInternalEdges()
        );

        for (label bndEdgei = 0; bndEdgei < nBoundaryEdges; ++bndEdgei)
        {
            edgeToBoundaryIndex.insert
            (
                patch().meshEdge(bndEdgei + nInternalEdges),
                edgeToBoundaryIndex.size()
            );

            // The attached patch face (there is only one):
            const label patchFacei = bndEdgeToPatchFace[bndEdgei][0];

            const label meshFacei = faceLabels_[patchFacei];

            const label patchId = pbm.whichPatch(meshFacei);

            bndEdgePatchPairs[bndEdgei].insert(patchId);
        }
    }

    // Pass 2:
    // - Add in first neighbour patch for the boundary edges
    // - examine all possible connecting neighbours

    findEdgePatchPairing
    (
        pbm,
        edgeToBoundaryIndex,
        bndEdgePatchPairs,
        edgeToBoundaryIndex.size()  // Number of places to still fill
    );


    // Nothing dangling if running in serial - can return already
    if (!Pstream::parRun())
    {
        return bndEdgePatchPairs;
    }

    // In parallel need to check for "dangling" edges, which are finiteArea
    // boundary edges that only exist on one side of a proc boundary.
    // Eg, proc boundary coincides with the outer extent of the finiteArea

    const globalMeshData& globalData = mesh().globalData();
    const indirectPrimitivePatch& cpp = globalData.coupledPatch();
    const mapDistribute& map = globalData.globalEdgeSlavesMap();

    // Construct coupled edge usage with all data
    List<unsigned char> coupledEdgesUsed(map.constructSize(), 0u);

    forAll(cpp.edges(), coupledEdgei)
    {
        const auto iter =
            edgeToBoundaryIndex.cfind(cpp.meshEdge(coupledEdgei));

        // Used from this side or other other side
        coupledEdgesUsed[coupledEdgei] = (iter.found() ? 1 : 2);
    }

    // Save the original (pre-sync) coupling state
    const List<unsigned char> coupledEdgesOrig(coupledEdgesUsed);

    globalData.syncData
    (
        coupledEdgesUsed,
        globalData.globalEdgeSlaves(),
        globalData.globalEdgeTransformedSlaves(),  // probably not used
        map,
        bitOrEqOp<unsigned char>()
    );


    // Check for one-sided edge coupling (coupled value == 3)
    // original == 1:
    // - coupled to a real finiteArea edge.
    // - receive a patch-pair value
    //
    // original == 2:
    // - a "dangled" edge. Information required for other procs.
    // - push a patch-pair value

    // Map edges (mesh numbering) back to a coupled index.
    // These are the edges to 'push' information for.

    EdgeMap<label> edgeToCoupledIndex;

    label nEdgesPull = 0;

    forAll(coupledEdgesUsed, coupledEdgei)
    {
        if (coupledEdgesUsed[coupledEdgei] == 3)
        {
            if (coupledEdgesOrig[coupledEdgei] == 1)
            {
                // Coupled side with finiteArea
                ++nEdgesPull;
            }
            else if (coupledEdgesOrig[coupledEdgei] == 2)
            {
                // Coupled side without finiteArea
                edgeToCoupledIndex.insert
                (
                    cpp.meshEdge(coupledEdgei),
                    coupledEdgei
                );
            }
        }
    }

    // Nothing to do - can return already
    if (returnReduce(edgeToCoupledIndex.empty(), andOp<bool>()))
    {
        return bndEdgePatchPairs;
    }

    // Data locations to pull
    labelList patchEdgeLabels(nEdgesPull);
    labelList coupledEdgeLabels(nEdgesPull);

    // Populate the locations
    {
        nEdgesPull = 0;

        forAll(cpp.edges(), coupledEdgei)
        {
            if
            (
                coupledEdgesUsed[coupledEdgei] == 3
             && coupledEdgesOrig[coupledEdgei] == 1
            )
            {
                // Pull this edge
                const auto iter =
                    edgeToBoundaryIndex.cfind(cpp.meshEdge(coupledEdgei));

                if (iter.found())
                {
                    patchEdgeLabels[nEdgesPull] = iter.val();
                    coupledEdgeLabels[nEdgesPull] = coupledEdgei;
                    ++nEdgesPull;
                }
                else
                {
                    // Should be impossible to fail here
                    FatalErrorInFunction
                        << "Failed on second lookup of "
                        << cpp.meshEdge(coupledEdgei) << nl
                        << abort(FatalError);
                }
            }
        }

        if (nEdgesPull != coupledEdgeLabels.size())
        {
            FatalErrorInFunction
                << "Failed lookup of some coupled edges" << nl
                << abort(FatalError);
        }
    }

    //- Construct edge sync with all data
    patchPairInfoList cppEdgeData(map.constructSize());

    // Fill in for 'push' locations. Only really interested in the owner
    // (corresponds to non-proc connection), but grab everything
    findEdgePatchPairing
    (
        pbm,
        edgeToCoupledIndex,
        cppEdgeData,
        2*edgeToCoupledIndex.size()  // Accept both sides?
    );


    const label nNonProcessor = pbm.nNonProcessor();

    // Adjust patch information to reflect dangling patch neighbour
    // Tag with -(value+2)
    forAllConstIters(edgeToCoupledIndex, iter)
    {
        const edge& e = iter.key();
        const label coupledEdgei = iter.val();

        patchPairInfo& pairing = cppEdgeData[coupledEdgei];
        const label ownerPatchId = pairing.first();

        // Some sanity checks
        if (ownerPatchId < 0)
        {
            FatalErrorInFunction
                << "Error finding dangling edge at "
                << cpp.points()[e.first()] << ' '
                << cpp.points()[e.second()] << nl
                << abort(FatalError);
        }
        else if (ownerPatchId >= nNonProcessor)
        {
            FatalErrorInFunction
                << "Cannot handle edge on processor-processor connection at "
                << cpp.points()[e.first()] << ' '
                << cpp.points()[e.second()] << nl
                << abort(FatalError);
        }

        combineDanglingEdge::setDangling(pairing, ownerPatchId);

        // TBD:
        // may wish to remember the corresponding proc number,
        // if we wish to bridge across 'fan-like' connections.
        //
        // pairing.setIndex(-(Pstream::myProcNo() + 2));
    }


    // Synchronize edge information

    const combineDanglingEdge edgeCombineOp(nNonProcessor);

    globalMeshData::syncData
    (
        cppEdgeData,
        globalData.globalEdgeSlaves(),
        globalData.globalEdgeTransformedSlaves(),  // probably not used
        map,
        edgeCombineOp
    );


    // Combine back from pushed cpp-edge data

    forAll(patchEdgeLabels, i)
    {
        patchPairInfo& pairing = bndEdgePatchPairs[patchEdgeLabels[i]];
        const patchPairInfo& other = cppEdgeData[coupledEdgeLabels[i]];

        edgeCombineOp(pairing, other);

        // Resolve special tagging
        combineDanglingEdge::correct(pairing);
    }

    return bndEdgePatchPairs;
}


Foam::PtrList<Foam::faPatch> Foam::faMesh::createPatchList
(
    const dictionary& bndDict,
    const word& emptyPatchName,
    const dictionary* defaultPatchDefinition
) const
{
    const polyBoundaryMesh& pbm = mesh().boundaryMesh();

    // Transcribe into patch definitions
    DynamicList<faPatchData> faPatchDefs(bndDict.size() + 4);
    for (const entry& dEntry : bndDict)
    {
        if (!dEntry.isDict())
        {
            WarningInFunction
                << "Not a dictionary entry: " << dEntry.name() << nl;
            continue;
        }
        const dictionary& patchDict = dEntry.dict();

        // Add entry
        faPatchDefs.append(faPatchData());

        auto& patchDef = faPatchDefs.last();
        patchDef.name_ = dEntry.keyword();
        patchDef.type_ = patchDict.get<word>("type");

        word patchName;

        // Optional: ownerPolyPatch
        if (patchDict.readIfPresent("ownerPolyPatch", patchName))
        {
            patchDef.ownerPolyPatchId_ = pbm.findPatchID(patchName);
            if (patchDef.ownerPolyPatchId_ < 0)
            {
                FatalErrorInFunction
                    << "ownerPolyPatch " << patchName << " not found"
                    << exit(FatalError);
            }
        }

        // Mandatory: neighbourPolyPatch
        patchDict.readEntry("neighbourPolyPatch", patchName);
        {
            patchDef.neighPolyPatchId_ = pbm.findPatchID(patchName);
            if (patchDef.neighPolyPatchId_ < 0)
            {
                FatalErrorInFunction
                    << "neighbourPolyPatch " << patchName << " not found"
                    << exit(FatalError);
            }
        }
    }

    // Additional empty placeholder patch?
    if (!emptyPatchName.empty())
    {
        faPatchDefs.append(faPatchData());

        auto& patchDef = faPatchDefs.last();
        patchDef.name_ = emptyPatchName;
        patchDef.type_ = "empty";
    }

    // Placeholder for any undefined edges
    const label undefPatchId = faPatchDefs.size();
    {
        faPatchDefs.append(faPatchData());

        auto& patchDef = faPatchDefs.last();
        patchDef.name_ = "undefined";
        patchDef.type_ = "patch";

        if (defaultPatchDefinition)
        {
            (*defaultPatchDefinition).readIfPresent("name", patchDef.name_);
            (*defaultPatchDefinition).readIfPresent("type", patchDef.type_);
        }
    }

    // ----------------------------------------------------------------------

    const label nInternalEdges = patch().nInternalEdges();
    const label nBoundaryEdges = patch().nBoundaryEdges();

    patchPairInfoList bndEdgePatchPairs(this->getBoundaryEdgePatchPairs());

    labelList bndEdgeFaPatchIDs(nBoundaryEdges, -1);

    for (label bndEdgei = 0; bndEdgei < nBoundaryEdges; ++bndEdgei)
    {
        const patchPairInfo& patchPair = bndEdgePatchPairs[bndEdgei];

        // Find first definition with a matching neighbour and
        // possibly with a matching owner.

        label bestPatchi = -1;

        for (label patchi = 0; patchi < faPatchDefs.size(); ++patchi)
        {
            const int match = faPatchDefs[patchi].matchPatchPair(patchPair);
            if (match == 3)
            {
                // Match (owner/neighbour) - done!
                bestPatchi = patchi;
                break;
            }
            else if (match == 2 && bestPatchi < 0)
            {
                // Match (neighbour) - keep looking for exact match
                bestPatchi = patchi;
            }
        }

        bndEdgeFaPatchIDs[bndEdgei] = bestPatchi;
    }


    // Extract which edges map to which patch
    // and set edgeLabels for each faPatch

    DynamicList<label> selectEdges(bndEdgeFaPatchIDs.size());

    for (label patchi = 0; patchi < faPatchDefs.size(); ++patchi)
    {
        auto& patchDef = faPatchDefs[patchi];

        selectEdges.clear();

        forAll(bndEdgeFaPatchIDs, bndEdgei)
        {
            if (bndEdgeFaPatchIDs[bndEdgei] == patchi)
            {
                selectEdges.append(bndEdgei + nInternalEdges);
            }
        }

        patchDef.edgeLabels_ = selectEdges;
    }

    // Check for undefined edges
    selectEdges.clear();

    forAll(bndEdgeFaPatchIDs, bndEdgei)
    {
        if (bndEdgeFaPatchIDs[bndEdgei] == -1)
        {
            selectEdges.append(bndEdgei + nInternalEdges);
        }
    }

    // Save the information
    faPatchDefs[undefPatchId].edgeLabels_ = selectEdges;

    bool hasUndefined = returnReduce(!selectEdges.empty(), orOp<bool>());

    if (hasUndefined)
    {
        // The initial edges to consider
        const labelList& undefinedEdges =
            faPatchDefs[undefPatchId].edgeLabels_;

        // Check for edges that butt against a processor (or other) patch

        labelList edgeNbrPolyPatch(undefinedEdges.size(), -1);
        forAll(edgeNbrPolyPatch, edgei)
        {
            const label patchEdgei = undefinedEdges[edgei];
            const label bndEdgei = (patchEdgei - nInternalEdges);

            edgeNbrPolyPatch[edgei] = bndEdgePatchPairs[bndEdgei].second();
        }

        // Categorize as processor/non-processor associations
        labelHashSet procPatchIDs;
        labelHashSet nonProcPatchIDs;

        for (const label polyPatchID : edgeNbrPolyPatch)
        {
            if (polyPatchID == -1)
            {
                nonProcPatchIDs.insert(polyPatchID);
            }
            else if
            (
                !nonProcPatchIDs.found(polyPatchID)
             && !procPatchIDs.found(polyPatchID)
            )
            {
                if (isA<processorPolyPatch>(pbm[polyPatchID]))
                {
                    procPatchIDs.insert(polyPatchID);
                }
                else
                {
                    nonProcPatchIDs.insert(polyPatchID);
                }
            }
        }

        // Select by processor association
        for (const label polyPatchID : procPatchIDs.sortedToc())
        {
            selectEdges.clear();

            forAll(edgeNbrPolyPatch, edgei)
            {
                if (edgeNbrPolyPatch[edgei] == polyPatchID)
                {
                    selectEdges.append(undefinedEdges[edgei]);
                }
            }

            faPatchDefs.append(faPatchData());

            auto& patchDef = faPatchDefs.last();
            patchDef.name_ = pbm[polyPatchID].name();
            patchDef.type_ = processorFaPatch::typeName;
            patchDef.neighPolyPatchId_ = polyPatchID;  // Needed for reorder

            const auto* ppp = isA<processorPolyPatch>(pbm[polyPatchID]);
            if (ppp)
            {
                patchDef.ownerProcId_ = ppp->myProcNo();
                patchDef.neighProcId_ = ppp->neighbProcNo();
            }

            patchDef.edgeLabels_ = selectEdges;
        }


        // Check for any remaining undefined edges
        selectEdges.clear();

        // Simply grab any/all (don't worry about which patch)
        if (!nonProcPatchIDs.empty())
        {
            forAll(edgeNbrPolyPatch, edgei)
            {
                const label polyPatchID = edgeNbrPolyPatch[edgei];
                if (nonProcPatchIDs.found(polyPatchID))
                {
                    selectEdges.append(undefinedEdges[edgei]);
                }
            }
        }

        // Complete the information
        faPatchDefs[undefPatchId].edgeLabels_ = selectEdges;

        hasUndefined = returnReduce(!selectEdges.empty(), orOp<bool>());
    }

    // Remove unnecessary entry
    if (!hasUndefined)
    {
        faPatchDefs.remove(undefPatchId);
    }

    for (auto& patchDef : faPatchDefs)
    {
        if (patchDef.coupled())
        {
            reorderProcEdges(patchDef, bndEdgePatchPairs);
            patchDef.neighPolyPatchId_ = -1; // No lookup of neighbour faces
        }
    }

    // Now convert list of definitions to list of patches

    label nPatches = 0;
    PtrList<faPatch> newPatches(faPatchDefs.size());

    for (faPatchData& patchDef : faPatchDefs)
    {
        newPatches.set
        (
            nPatches,
            faPatch::New
            (
                patchDef.name(),        // name
                patchDef.dict(false),   // withEdgeLabels == false
                nPatches,               // index
                boundary()
            )
        );

        // Transfer edge labels
        newPatches[nPatches].resetEdges(std::move(patchDef.edgeLabels_));
        ++nPatches;
    }

    return newPatches;
}


// ************************************************************************* //
