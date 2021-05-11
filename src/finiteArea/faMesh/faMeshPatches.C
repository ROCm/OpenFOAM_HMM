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
#include "edgeHashes.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faMesh::reorderProcEdges
(
    faPatchData& patchDef,
    const labelUList& meshEdges
) const
{
    if (!patchDef.coupled() || patchDef.edgeLabels_.empty())
    {
        return;
    }

    const polyBoundaryMesh& pbm = mesh().boundaryMesh();
    const labelListList& edgeFaces = mesh().edgeFaces();

    // Reorder processor edges using order of neighbour processorPolyPatch

    const label procPatchID = patchDef.neighPolyPatchId_;
    const label nProcEdges = patchDef.edgeLabels_.size();

    labelList procFaces(nProcEdges, -1);

    forAll(procFaces, edgei)
    {
        const label localEdgei = patchDef.edgeLabels_[edgei];
        const label meshEdgei = meshEdges[localEdgei];

        for (const label meshFacei : edgeFaces[meshEdgei])
        {
            if
            (
                !faceLabels_.found(meshFacei)
             && (procPatchID == pbm.whichPatch(meshFacei))
            )
            {
                // The edge's proc-face
                procFaces[edgei] = meshFacei;
                break;
            }
        }
    }

    // Ascending proc-face numbering
    const labelList sortIndices(Foam::sortedOrder(procFaces));

    const labelList& oldEdgeLabels = patchDef.edgeLabels_;
    labelList newEdgeLabels(oldEdgeLabels.size());

    // Most of the time, an individual proc-face will only be singly
    // attached to the finite-area patch. In rarer case, there could
    // multiple connections. For these cases, need to walk the face
    // edges - the direction depends on owner vs neighbour side.

    EdgeMap<label> multihit;

    for (label edgei = 0; edgei < nProcEdges; /*nil*/)
    {
        const label procFacei = procFaces[sortIndices[edgei]];

        // Find all identical faces
        label endEdgei = edgei + 1;  // one beyond
        while
        (
            (endEdgei < nProcEdges)
         && (procFacei == procFaces[sortIndices[endEdgei]])
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

            // Map from global edge to local edgeId
            for (label i = edgei; i < endEdgei; ++i)
            {
                label localEdgei = oldEdgeLabels[sortIndices[i]];
                label meshEdgei = meshEdges[localEdgei];

                multihit.insert
                (
                    mesh().edges()[meshEdgei],
                    localEdgei
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

            const face& f = mesh().faces()[procFacei];

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
                    << " edges for face: " << procFacei
                    << " ... indicates serious geometry issue" << nl
                    << multihit << nl
                    << abort(FatalError);
            }
            if (!multihit.empty())
            {
                FatalErrorInFunction
                    << "Missed edges for face: " << procFacei
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
        "",             // Name for empty patch placeholder
        &onePatchDict   // Definitions for defaultPatch
    );
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

        const word ownName(patchDict.get<word>("ownerPolyPatch"));
        const word neiName(patchDict.get<word>("neighbourPolyPatch"));

        patchDef.ownerPolyPatchId_ = pbm.findPatchID(ownName);
        patchDef.neighPolyPatchId_ = pbm.findPatchID(neiName);

        if (patchDef.ownerPolyPatchId_ < 0)
        {
            FatalErrorInFunction
                << "ownerPolyPatch " << ownName << " not found"
                << exit(FatalError);
        }
        if (patchDef.neighPolyPatchId_ < 0)
        {
            FatalErrorInFunction
                << "neighbourPolyPatch " << neiName << " not found"
                << exit(FatalError);
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

    // Determine faPatch ID for each boundary edge.
    // Result is in the bndEdgeFaPatchIDs list

    const labelList meshEdges
    (
        patch().meshEdges(mesh().edges(), mesh().pointEdges())
    );

    const labelListList& edgeFaces = mesh().edgeFaces();

    const label nInternalEdges = patch().nInternalEdges();
    const label nBoundaryEdges = patch().nBoundaryEdges();

    labelList bndEdgeFaPatchIDs(nBoundaryEdges, -1);

    for (label bndEdgei = 0; bndEdgei < nBoundaryEdges; ++bndEdgei)
    {
        const label patchEdgei = meshEdges[bndEdgei + nInternalEdges];

        // Use 'edge' for accounting
        edge curEdgePatchPair;

        for (const label meshFacei : edgeFaces[patchEdgei])
        {
            const label polyPatchID = pbm.whichPatch(meshFacei);

            if (polyPatchID != -1)
            {
                curEdgePatchPair.insert(polyPatchID);
            }
        }


        if (curEdgePatchPair.valid())
        {
            // Non-negative, unique pairing
            // - find corresponding definition

            for (label patchi = 0; patchi < faPatchDefs.size(); ++patchi)
            {
                if (faPatchDefs[patchi].foundPatchPair(curEdgePatchPair))
                {
                    bndEdgeFaPatchIDs[bndEdgei] = patchi;
                    break;
                }
            }
        }
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
            const label localEdgei = undefinedEdges[edgei];
            const label meshEdgei = meshEdges[localEdgei];
            label polyPatchID;

            for (const label meshFacei : edgeFaces[meshEdgei])
            {
                if
                (
                    !faceLabels_.found(meshFacei)
                 && (polyPatchID = pbm.whichPatch(meshFacei)) != -1
                )
                {
                    // Found the edge's off-patch neighbour face
                    edgeNbrPolyPatch[edgei] = polyPatchID;
                    break;
                }
            }
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
            reorderProcEdges(patchDef, meshEdges);
            patchDef.neighPolyPatchId_ = -1; // No longer required + confusing
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
