/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "processorPolyPatch.H"
#include "PstreamBuffers.H"
#include "lduSchedule.H"
#include "globalMeshData.H"
#include "stringListOps.H"
#include "DynamicList.H"
#include "PtrListOps.H"
#include "edgeHashes.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyBoundaryMesh, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::polyBoundaryMesh::hasGroupIDs() const
{
    if (groupIDsPtr_)
    {
        // Use existing cache
        return !groupIDsPtr_->empty();
    }

    const polyPatchList& patches = *this;

    for (const polyPatch& p : patches)
    {
        if (!p.inGroups().empty())
        {
            return true;
        }
    }

    return false;
}


void Foam::polyBoundaryMesh::calcGroupIDs() const
{
    if (groupIDsPtr_)
    {
        return;  // Or FatalError
    }

    groupIDsPtr_.reset(new HashTable<labelList>(16));
    auto& groupLookup = *groupIDsPtr_;

    const polyPatchList& patches = *this;

    forAll(patches, patchi)
    {
        const wordList& groups = patches[patchi].inGroups();

        for (const word& groupName : groups)
        {
            groupLookup(groupName).append(patchi);
        }
    }

    // Remove groups that clash with patch names
    forAll(patches, patchi)
    {
        if (groupLookup.erase(patches[patchi].name()))
        {
            WarningInFunction
                << "Removed group '" << patches[patchi].name()
                << "' which clashes with patch " << patchi
                << " of the same name."
                << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyBoundaryMesh::polyBoundaryMesh
(
    const IOobject& io,
    const polyMesh& mesh
)
:
    polyPatchList(),
    regIOobject(io),
    mesh_(mesh)
{
    if
    (
        readOpt() == IOobject::MUST_READ
     || readOpt() == IOobject::MUST_READ_IF_MODIFIED
    )
    {
        // Warn for MUST_READ_IF_MODIFIED
        warnNoRereading<polyBoundaryMesh>();

        polyPatchList& patches = *this;

        // Read polyPatchList
        Istream& is = readStream(typeName);

        PtrList<entry> patchEntries(is);
        patches.setSize(patchEntries.size());

        forAll(patches, patchi)
        {
            patches.set
            (
                patchi,
                polyPatch::New
                (
                    patchEntries[patchi].keyword(),
                    patchEntries[patchi].dict(),
                    patchi,
                    *this
                )
            );
        }

        is.check(FUNCTION_NAME);

        close();
    }
}


Foam::polyBoundaryMesh::polyBoundaryMesh
(
    const IOobject& io,
    const polyMesh& pm,
    const label size
)
:
    polyPatchList(size),
    regIOobject(io),
    mesh_(pm)
{}


Foam::polyBoundaryMesh::polyBoundaryMesh
(
    const IOobject& io,
    const polyMesh& pm,
    const polyPatchList& ppl
)
:
    polyPatchList(),
    regIOobject(io),
    mesh_(pm)
{
    if
    (
        (this->readOpt() == IOobject::READ_IF_PRESENT && this->headerOk())
     || this->readOpt() == IOobject::MUST_READ
     || this->readOpt() == IOobject::MUST_READ_IF_MODIFIED
    )
    {
        // Warn for MUST_READ_IF_MODIFIED
        warnNoRereading<polyBoundaryMesh>();

        polyPatchList& patches = *this;

        // Read polyPatchList
        Istream& is = readStream(typeName);

        PtrList<entry> patchEntries(is);
        patches.setSize(patchEntries.size());

        forAll(patches, patchi)
        {
            patches.set
            (
                patchi,
                polyPatch::New
                (
                    patchEntries[patchi].keyword(),
                    patchEntries[patchi].dict(),
                    patchi,
                    *this
                )
            );
        }

        is.check(FUNCTION_NAME);

        close();
    }
    else
    {
        polyPatchList& patches = *this;
        patches.setSize(ppl.size());
        forAll(patches, patchi)
        {
            patches.set(patchi, ppl[patchi].clone(*this).ptr());
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

void Foam::polyBoundaryMesh::clearGeom()
{
    polyPatchList& patches = *this;

    for (polyPatch& p : patches)
    {
        p.clearGeom();
    }
}


void Foam::polyBoundaryMesh::clearAddressing()
{
    neighbourEdgesPtr_.clear();
    patchIDPtr_.clear();
    groupIDsPtr_.clear();

    polyPatchList& patches = *this;

    for (polyPatch& p : patches)
    {
        p.clearAddressing();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyBoundaryMesh::calcGeometry()
{
    PstreamBuffers pBufs(Pstream::defaultCommsType);

    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        forAll(*this, patchi)
        {
            operator[](patchi).initGeometry(pBufs);
        }

        pBufs.finishedSends();

        forAll(*this, patchi)
        {
            operator[](patchi).calcGeometry(pBufs);
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = mesh().globalData().patchSchedule();

        // Dummy.
        pBufs.finishedSends();

        for (const auto& patchEval : patchSchedule)
        {
            const label patchi = patchEval.patch;

            if (patchEval.init)
            {
                operator[](patchi).initGeometry(pBufs);
            }
            else
            {
                operator[](patchi).calcGeometry(pBufs);
            }
        }
    }
}


Foam::UPtrList<const Foam::labelUList>
Foam::polyBoundaryMesh::faceCells() const
{
    const polyPatchList& patches = *this;

    UPtrList<const labelUList> list(patches.size());

    forAll(patches, patchi)
    {
        list.set(patchi, &patches[patchi].faceCells());
    }

    return list;
}


const Foam::List<Foam::labelPairList>&
Foam::polyBoundaryMesh::neighbourEdges() const
{
    if (Pstream::parRun())
    {
        WarningInFunction
            << "Neighbour edge addressing not correct across parallel"
            << " boundaries." << endl;
    }

    if (!neighbourEdgesPtr_)
    {
        neighbourEdgesPtr_.reset(new List<labelPairList>(size()));
        auto& neighbourEdges = *neighbourEdgesPtr_;

        // Initialize.
        label nEdgePairs = 0;
        forAll(*this, patchi)
        {
            const polyPatch& pp = operator[](patchi);

            neighbourEdges[patchi].setSize(pp.nEdges() - pp.nInternalEdges());

            for (labelPair& edgeInfo : neighbourEdges[patchi])
            {
                edgeInfo[0] = -1;
                edgeInfo[1] = -1;
            }

            nEdgePairs += pp.nEdges() - pp.nInternalEdges();
        }

        // From mesh edge (expressed as a point pair so as not to construct
        // point addressing) to patch + relative edge index.
        EdgeMap<labelPair> pointsToEdge(nEdgePairs);

        forAll(*this, patchi)
        {
            const polyPatch& pp = operator[](patchi);

            const edgeList& edges = pp.edges();

            for
            (
                label edgei = pp.nInternalEdges();
                edgei < edges.size();
                edgei++
            )
            {
                // Edge in patch local points
                const edge& e = edges[edgei];

                // Edge in mesh points.
                edge meshEdge(pp.meshPoints()[e[0]], pp.meshPoints()[e[1]]);

                auto fnd = pointsToEdge.find(meshEdge);

                if (!fnd.found())
                {
                    // First occurrence of mesh edge. Store patch and my
                    // local index.
                    pointsToEdge.insert
                    (
                        meshEdge,
                        labelPair
                        (
                            patchi,
                            edgei - pp.nInternalEdges()
                        )
                    );
                }
                else
                {
                    // Second occurrence. Store.
                    const labelPair& edgeInfo = fnd.val();

                    neighbourEdges[patchi][edgei - pp.nInternalEdges()] =
                        edgeInfo;

                    neighbourEdges[edgeInfo[0]][edgeInfo[1]]
                         = labelPair(patchi, edgei - pp.nInternalEdges());

                    // Found all two occurrences of this edge so remove from
                    // hash to save space. Note that this will give lots of
                    // problems if the polyBoundaryMesh is multiply connected.
                    pointsToEdge.erase(meshEdge);
                }
            }
        }

        if (pointsToEdge.size())
        {
            FatalErrorInFunction
                << "Not all boundary edges of patches match up." << nl
                << "Is the outside of your mesh multiply connected?"
                << abort(FatalError);
        }

        forAll(*this, patchi)
        {
            const polyPatch& pp = operator[](patchi);

            const labelPairList& nbrEdges = neighbourEdges[patchi];

            forAll(nbrEdges, i)
            {
                const labelPair& edgeInfo = nbrEdges[i];

                if (edgeInfo[0] == -1 || edgeInfo[1] == -1)
                {
                    const label edgei = pp.nInternalEdges() + i;
                    const edge& e = pp.edges()[edgei];

                    FatalErrorInFunction
                        << "Not all boundary edges of patches match up." << nl
                        << "Edge " << edgei << " on patch " << pp.name()
                        << " end points " << pp.localPoints()[e[0]] << ' '
                        << pp.localPoints()[e[1]] << " is not matched to an"
                        << " edge on any other patch." << nl
                        << "Is the outside of your mesh multiply connected?"
                        << abort(FatalError);
                }
            }
        }
    }

    return *neighbourEdgesPtr_;
}


const Foam::labelList& Foam::polyBoundaryMesh::patchID() const
{
    if (!patchIDPtr_)
    {
        patchIDPtr_.reset(new labelList(mesh_.nBoundaryFaces()));
        labelList& list = *patchIDPtr_;

        const polyPatchList& patches = *this;

        forAll(patches, patchi)
        {
            SubList<label>
            (
                list,
                patches[patchi].size(),
                (patches[patchi].start() - mesh_.nInternalFaces())
            ) = patchi;
        }
    }

    return *patchIDPtr_;
}


const Foam::HashTable<Foam::labelList>&
Foam::polyBoundaryMesh::groupPatchIDs() const
{
    if (!groupIDsPtr_)
    {
        calcGroupIDs();
    }

    return *groupIDsPtr_;
}


void Foam::polyBoundaryMesh::setGroup
(
    const word& groupName,
    const labelUList& patchIDs
)
{
    groupIDsPtr_.clear();

    polyPatchList& patches = *this;

    boolList donePatch(patches.size(), false);

    // Add to specified patches
    for (const label patchi : patchIDs)
    {
        patches[patchi].inGroups().appendUniq(groupName);
        donePatch[patchi] = true;
    }

    // Remove from other patches
    forAll(patches, patchi)
    {
        if (!donePatch[patchi])
        {
            wordList& groups = patches[patchi].inGroups();

            if (groups.found(groupName))
            {
                label newi = 0;
                forAll(groups, i)
                {
                    if (groups[i] != groupName)
                    {
                        groups[newi++] = groups[i];
                    }
                }
                groups.resize(newi);
            }
        }
    }
}


Foam::label Foam::polyBoundaryMesh::nNonProcessor() const
{
    const polyPatchList& patches = *this;

    label nonProc = 0;

    for (const polyPatch& p : patches)
    {
        if (isA<processorPolyPatch>(p))
        {
            break;
        }

        ++nonProc;
    }

    return nonProc;
}


Foam::wordList Foam::polyBoundaryMesh::names() const
{
    return PtrListOps::get<word>(*this, nameOp<polyPatch>());
}


Foam::wordList Foam::polyBoundaryMesh::types() const
{
    return PtrListOps::get<word>(*this, typeOp<polyPatch>());
}


Foam::wordList Foam::polyBoundaryMesh::physicalTypes() const
{
    return
        PtrListOps::get<word>
        (
            *this,
            [](const polyPatch& p) { return p.physicalType(); }
        );
}


Foam::labelList Foam::polyBoundaryMesh::patchStarts() const
{
    return
        PtrListOps::get<label>
        (
            *this,
            [](const polyPatch& p) { return p.start(); }
        );
}


Foam::labelList Foam::polyBoundaryMesh::patchSizes() const
{
    return
        PtrListOps::get<label>
        (
            *this,
            [](const polyPatch& p) { return p.size(); }
        );
}


Foam::List<Foam::labelRange> Foam::polyBoundaryMesh::patchRanges() const
{
    return
        PtrListOps::get<labelRange>
        (
            *this,
            [](const polyPatch& p) { return p.range(); }
        );
}


Foam::label Foam::polyBoundaryMesh::start() const
{
    return mesh_.nInternalFaces();
}


Foam::label Foam::polyBoundaryMesh::nFaces() const
{
    return mesh_.nBoundaryFaces();
}


Foam::labelRange Foam::polyBoundaryMesh::range() const
{
    return labelRange(mesh_.nInternalFaces(), mesh_.nBoundaryFaces());
}


Foam::labelRange Foam::polyBoundaryMesh::range(const label patchi) const
{
    if (patchi < 0)
    {
        return labelRange(mesh_.nInternalFaces(), 0);
    }

    // Will fail if patchi >= size()
    return (*this)[patchi].range();
}


Foam::labelList Foam::polyBoundaryMesh::indices
(
    const wordRe& matcher,
    const bool useGroups
) const
{
    if (matcher.empty())
    {
        return labelList();
    }

    // Only check groups if requested and they exist
    const bool checkGroups = (useGroups && this->hasGroupIDs());

    labelHashSet ids;

    if (matcher.isPattern())
    {
        if (checkGroups)
        {
            const auto& groupLookup = groupPatchIDs();
            forAllConstIters(groupLookup, iter)
            {
                if (matcher.match(iter.key()))
                {
                    // Hash ids associated with the group
                    ids.insert(iter.val());
                }
            }
        }

        if (ids.empty())
        {
            return PtrListOps::findMatching(*this, matcher);
        }
        else
        {
            ids.insert(PtrListOps::findMatching(*this, matcher));
        }
    }
    else
    {
        // Literal string.
        // Special version of above for reduced memory footprint.

        const label patchId = PtrListOps::firstMatching(*this, matcher);

        if (patchId >= 0)
        {
            return labelList(one{}, patchId);
        }
        else if (checkGroups)
        {
            const auto iter = groupPatchIDs().cfind(matcher);

            if (iter.found())
            {
                // Hash ids associated with the group
                ids.insert(iter.val());
            }
        }
    }

    return ids.sortedToc();
}


Foam::labelList Foam::polyBoundaryMesh::indices
(
    const wordRes& matcher,
    const bool useGroups
) const
{
    if (matcher.empty())
    {
        return labelList();
    }
    else if (matcher.size() == 1)
    {
        return this->indices(matcher.first(), useGroups);
    }

    labelHashSet ids;

    // Only check groups if requested and they exist
    if (useGroups && this->hasGroupIDs())
    {
        ids.resize(2*this->size());

        const auto& groupLookup = groupPatchIDs();
        forAllConstIters(groupLookup, iter)
        {
            if (matcher.match(iter.key()))
            {
                // Hash ids associated with the group
                ids.insert(iter.val());
            }
        }
    }

    if (ids.empty())
    {
        return PtrListOps::findMatching(*this, matcher);
    }
    else
    {
        ids.insert(PtrListOps::findMatching(*this, matcher));
    }

    return ids.sortedToc();
}


Foam::label Foam::polyBoundaryMesh::findIndex(const wordRe& key) const
{
    if (key.empty())
    {
        return -1;
    }
    return PtrListOps::firstMatching(*this, key);
}


Foam::label Foam::polyBoundaryMesh::findPatchID
(
    const word& patchName,
    bool allowNotFound
) const
{
    if (patchName.empty())
    {
        return -1;
    }

    const label patchId = PtrListOps::firstMatching(*this, patchName);

    if (patchId >= 0)
    {
        return patchId;
    }

    if (!allowNotFound)
    {
        FatalErrorInFunction
            << "Patch '" << patchName << "' not found. "
            << "Available patch names";

        if (polyMesh::defaultRegion != mesh_.name())
        {
            FatalError
                << " in region '" << mesh_.name() << "'";
        }

        FatalError
            << " include: " << names() << endl
            << exit(FatalError);
    }

    // Patch not found
    if (debug)
    {
        Pout<< "label polyBoundaryMesh::findPatchID(const word&) const"
            << "Patch named " << patchName << " not found.  "
            << "Available patch names: " << names() << endl;
    }

    // Not found, return -1
    return -1;
}


Foam::label Foam::polyBoundaryMesh::whichPatch(const label faceIndex) const
{
    // Find out which patch the current face belongs to by comparing label
    // with patch start labels.
    // If the face is internal, return -1;
    // if it is off the end of the list, abort
    if (faceIndex < mesh().nInternalFaces())
    {
        return -1;
    }
    else if (faceIndex >= mesh().nFaces())
    {
        FatalErrorInFunction
            << "Face " << faceIndex
            << " out of bounds. Number of geometric faces " << mesh().nFaces()
            << abort(FatalError);
    }


    // Patches are ordered, use binary search

    const polyPatchList& patches = *this;

    const label patchi =
        findLower
        (
            patches,
            faceIndex,
            0,
            // Must include the start in the comparison
            [](const polyPatch& p, label val) { return (p.start() <= val); }
        );

    if (patchi < 0 || !patches[patchi].range().found(faceIndex))
    {
        // If not in any of above, it is trouble!
        FatalErrorInFunction
            << "Face " << faceIndex << " not found in any of the patches "
            << flatOutput(names()) << nl
            << "The patches appear to be inconsistent with the mesh :"
            << " internalFaces:" << mesh().nInternalFaces()
            << " total number of faces:" << mesh().nFaces()
            << abort(FatalError);

        return -1;
    }

    return patchi;
}


Foam::labelHashSet Foam::polyBoundaryMesh::patchSet
(
    const UList<wordRe>& patchNames,
    const bool warnNotFound,
    const bool useGroups
) const
{
    const wordList allPatchNames(this->names());
    labelHashSet ids(2*this->size());

    // Only check groups if requested and they exist
    const bool checkGroups = (useGroups && this->hasGroupIDs());

    for (const wordRe& matcher : patchNames)
    {
        labelList matchIndices = findMatchingStrings(matcher, allPatchNames);
        ids.insert(matchIndices);

        bool missed = matchIndices.empty();

        if (missed && checkGroups)
        {
            // Check group names
            if (matcher.isPattern())
            {
                forAllConstIters(groupPatchIDs(), iter)
                {
                    if (matcher.match(iter.key()))
                    {
                        // Hash ids associated with the group
                        ids.insert(iter.val());
                        missed = false;
                    }
                }
            }
            else
            {
                const auto iter = groupPatchIDs().cfind(matcher);

                if (iter.found())
                {
                    // Hash ids associated with the group
                    ids.insert(iter.val());
                    missed = false;
                }
            }
        }

        if (missed && warnNotFound)
        {
            if (checkGroups)
            {
                WarningInFunction
                    << "Cannot find any patch or group names matching "
                    << matcher << endl;
            }
            else
            {
                WarningInFunction
                    << "Cannot find any patch names matching "
                    << matcher << endl;
            }
        }
    }

    return ids;
}


void Foam::polyBoundaryMesh::matchGroups
(
    const labelUList& patchIDs,
    wordList& groups,
    labelHashSet& nonGroupPatches
) const
{
    // Current matched groups
    DynamicList<word> matchedGroups(1);

    // Current set of unmatched patches
    nonGroupPatches = labelHashSet(patchIDs);

    const HashTable<labelList>& groupLookup = this->groupPatchIDs();
    forAllConstIters(groupLookup, iter)
    {
        // Store currently unmatched patches so we can restore
        labelHashSet oldNonGroupPatches(nonGroupPatches);

        // Match by deleting patches in group from the current set and seeing
        // if all have been deleted.
        labelHashSet groupPatchSet(iter.val());

        label nMatch = nonGroupPatches.erase(groupPatchSet);

        if (nMatch == groupPatchSet.size())
        {
            matchedGroups.append(iter.key());
        }
        else if (nMatch != 0)
        {
            // No full match. Undo.
            nonGroupPatches.transfer(oldNonGroupPatches);
        }
    }

    groups.transfer(matchedGroups);
}


bool Foam::polyBoundaryMesh::checkParallelSync(const bool report) const
{
    if (!Pstream::parRun())
    {
        return false;
    }

    const polyBoundaryMesh& bm = *this;

    bool hasError = false;

    // Collect non-proc patches and check proc patches are last.
    wordList names(bm.size());
    wordList types(bm.size());

    label nonProci = 0;

    forAll(bm, patchi)
    {
        if (!isA<processorPolyPatch>(bm[patchi]))
        {
            if (nonProci != patchi)
            {
                // There is processor patch in between normal patches.
                hasError = true;

                if (debug || report)
                {
                    Pout<< " ***Problem with boundary patch " << patchi
                        << " named " << bm[patchi].name()
                        << " of type " <<  bm[patchi].type()
                        << ". The patch seems to be preceeded by processor"
                        << " patches. This is can give problems."
                        << endl;
                }
            }
            else
            {
                names[nonProci] = bm[patchi].name();
                types[nonProci] = bm[patchi].type();
                nonProci++;
            }
        }
    }
    names.setSize(nonProci);
    types.setSize(nonProci);

    List<wordList> allNames(Pstream::nProcs());
    allNames[Pstream::myProcNo()] = names;
    Pstream::gatherList(allNames);
    Pstream::scatterList(allNames);

    List<wordList> allTypes(Pstream::nProcs());
    allTypes[Pstream::myProcNo()] = types;
    Pstream::gatherList(allTypes);
    Pstream::scatterList(allTypes);

    // Have every processor check but print error on master
    // (in processor sequence).

    for (const int proci : Pstream::subProcs())
    {
        if
        (
            (allNames[proci] != allNames.first())
         || (allTypes[proci] != allTypes.first())
        )
        {
            hasError = true;

            if (debug || (report && Pstream::master()))
            {
                Info<< " ***Inconsistent patches across processors, "
                       "processor 0 has patch names:"
                    << allNames.first()
                    << " patch types:" << allTypes.first()
                    << " processor " << proci
                    << " has patch names:" << allNames[proci]
                    << " patch types:" << allTypes[proci]
                    << endl;
            }
        }
    }

    return hasError;
}


bool Foam::polyBoundaryMesh::checkDefinition(const bool report) const
{
    label nextPatchStart = mesh().nInternalFaces();
    const polyBoundaryMesh& bm = *this;

    bool hasError = false;

    wordHashSet patchNames(2*this->size());

    forAll(bm, patchi)
    {
        if (bm[patchi].start() != nextPatchStart && !hasError)
        {
            hasError = true;

            Info<< " ****Problem with boundary patch " << patchi
                << " named " << bm[patchi].name()
                << " of type " <<  bm[patchi].type()
                << ". The patch should start on face no " << nextPatchStart
                << " and the patch specifies " << bm[patchi].start()
                << "." << endl
                << "Possibly consecutive patches have this same problem."
                << " Suppressing future warnings." << endl;
        }

        if (!patchNames.insert(bm[patchi].name()) && !hasError)
        {
            hasError = true;

            Info<< " ****Duplicate boundary patch " << patchi
                << " named " << bm[patchi].name()
                << " of type " <<  bm[patchi].type()
                << "." << endl
                << "Suppressing future warnings." << endl;
        }

        nextPatchStart += bm[patchi].size();
    }

    reduce(hasError, orOp<bool>());

    if (debug || report)
    {
        if (hasError)
        {
            Pout<< " ***Boundary definition is in error." << endl;
        }
        else
        {
            Info<< "    Boundary definition OK." << endl;
        }
    }

    return hasError;
}


void Foam::polyBoundaryMesh::movePoints(const pointField& p)
{
    PstreamBuffers pBufs(Pstream::defaultCommsType);

    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        forAll(*this, patchi)
        {
            operator[](patchi).initMovePoints(pBufs, p);
        }

        pBufs.finishedSends();

        forAll(*this, patchi)
        {
            operator[](patchi).movePoints(pBufs, p);
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = mesh().globalData().patchSchedule();

        // Dummy.
        pBufs.finishedSends();

        for (const auto& patchEval : patchSchedule)
        {
            const label patchi = patchEval.patch;

            if (patchEval.init)
            {
                operator[](patchi).initMovePoints(pBufs, p);
            }
            else
            {
                operator[](patchi).movePoints(pBufs, p);
            }
        }
    }
}


void Foam::polyBoundaryMesh::updateMesh()
{
    neighbourEdgesPtr_.clear();
    patchIDPtr_.clear();
    groupIDsPtr_.clear();

    PstreamBuffers pBufs(Pstream::defaultCommsType);

    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        forAll(*this, patchi)
        {
            operator[](patchi).initUpdateMesh(pBufs);
        }

        pBufs.finishedSends();

        forAll(*this, patchi)
        {
            operator[](patchi).updateMesh(pBufs);
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = mesh().globalData().patchSchedule();

        // Dummy.
        pBufs.finishedSends();

        for (const auto& patchEval : patchSchedule)
        {
            const label patchi = patchEval.patch;

            if (patchEval.init)
            {
                operator[](patchi).initUpdateMesh(pBufs);
            }
            else
            {
                operator[](patchi).updateMesh(pBufs);
            }
        }
    }
}


void Foam::polyBoundaryMesh::reorder
(
    const labelUList& oldToNew,
    const bool validBoundary
)
{
    // Change order of patches
    polyPatchList::reorder(oldToNew);

    // Adapt indices
    polyPatchList& patches = *this;

    forAll(patches, patchi)
    {
        patches[patchi].index() = patchi;
    }

    if (validBoundary)
    {
        updateMesh();
    }
}


bool Foam::polyBoundaryMesh::writeData(Ostream& os) const
{
    const polyPatchList& patches = *this;

    os  << patches.size() << nl << token::BEGIN_LIST << incrIndent << nl;

    for (const polyPatch& pp : patches)
    {
        os.beginBlock(pp.name());
        os  << pp;
        os.endBlock();
    }

    os  << decrIndent << token::END_LIST;

    os.check(FUNCTION_NAME);
    return os.good();
}


bool Foam::polyBoundaryMesh::writeObject
(
    IOstreamOption streamOpt,
    const bool valid
) const
{
    streamOpt.compression(IOstreamOption::UNCOMPRESSED);
    return regIOobject::writeObject(streamOpt, valid);
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

const Foam::polyPatch& Foam::polyBoundaryMesh::operator[]
(
    const word& patchName
) const
{
    const label patchi = findPatchID(patchName);

    if (patchi < 0)
    {
        FatalErrorInFunction
            << "Patch named " << patchName << " not found." << nl
            << "Available patch names: " << names() << endl
            << abort(FatalError);
    }

    return operator[](patchi);
}


Foam::polyPatch& Foam::polyBoundaryMesh::operator[]
(
    const word& patchName
)
{
    const label patchi = findPatchID(patchName);

    if (patchi < 0)
    {
        FatalErrorInFunction
            << "Patch named " << patchName << " not found." << nl
            << "Available patch names: " << names() << endl
            << abort(FatalError);
    }

    return operator[](patchi);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const polyBoundaryMesh& pbm)
{
    pbm.writeData(os);
    return os;
}


// ************************************************************************* //
