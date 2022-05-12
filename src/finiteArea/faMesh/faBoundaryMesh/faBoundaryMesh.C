/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "faBoundaryMesh.H"
#include "faMesh.H"
#include "globalIndex.H"
#include "primitiveMesh.H"
#include "processorFaPatch.H"
#include "PtrListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faBoundaryMesh, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::faBoundaryMesh::hasGroupIDs() const
{
    if (groupIDsPtr_)
    {
        // Use existing cache
        return !groupIDsPtr_->empty();
    }

    const faPatchList& patches = *this;

    for (const faPatch& p : patches)
    {
        if (!p.inGroups().empty())
        {
            return true;
        }
    }

    return false;
}


void Foam::faBoundaryMesh::calcGroupIDs() const
{
    if (groupIDsPtr_)
    {
        return;  // Or FatalError
    }

    groupIDsPtr_.reset(new HashTable<labelList>(16));
    auto& groupLookup = *groupIDsPtr_;

    const faPatchList& patches = *this;

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


bool Foam::faBoundaryMesh::readContents(const bool allowReadIfPresent)
{
    if
    (
        readOpt() == IOobject::MUST_READ
     || readOpt() == IOobject::MUST_READ_IF_MODIFIED
     ||
        (
            allowReadIfPresent
         && (readOpt() == IOobject::READ_IF_PRESENT && headerOk())
        )
    )
    {
        // Warn for MUST_READ_IF_MODIFIED
        warnNoRereading<faBoundaryMesh>();

        faPatchList& patches = *this;

        // Read faPatch list
        Istream& is = readStream(typeName);

        // Read patches as entries
        PtrList<entry> patchEntries(is);
        patches.resize(patchEntries.size());

        // Transcribe
        forAll(patches, patchi)
        {
            patches.set
            (
                patchi,
                faPatch::New
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
        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faBoundaryMesh::faBoundaryMesh
(
    const IOobject& io,
    const faMesh& mesh
)
:
    faPatchList(),
    regIOobject(io),
    mesh_(mesh)
{
    readContents(false);  // READ_IF_PRESENT allowed: False
}


Foam::faBoundaryMesh::faBoundaryMesh
(
    const IOobject& io,
    const faMesh& pm,
    const label size
)
:
    faPatchList(size),
    regIOobject(io),
    mesh_(pm)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faBoundaryMesh::calcGeometry()
{
    // processorFaPatch geometry triggers calculation of pointNormals.
    // This uses parallel comms and hence will not be trigggered
    // on processors that do not have a processorFaPatch so instead
    // force construction.
    (void)mesh_.pointAreaNormals();

    PstreamBuffers pBufs(Pstream::defaultCommsType);

    if
    (
        pBufs.commsType() == Pstream::commsTypes::blocking
     || pBufs.commsType() == Pstream::commsTypes::nonBlocking
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
    else if (pBufs.commsType() == Pstream::commsTypes::scheduled)
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
Foam::faBoundaryMesh::edgeLabels() const
{
    const faPatchList& patches = *this;

    UPtrList<const labelUList> list(patches.size());

    forAll(list, patchi)
    {
        list.set(patchi, &patches[patchi].edgeLabels());
    }

    return list;
}


Foam::UPtrList<const Foam::labelUList>
Foam::faBoundaryMesh::edgeFaces() const
{
    const faPatchList& patches = *this;

    UPtrList<const labelUList> list(patches.size());

    forAll(list, patchi)
    {
        list.set(patchi, &patches[patchi].edgeFaces());
    }

    return list;
}


Foam::lduInterfacePtrsList Foam::faBoundaryMesh::interfaces() const
{
    const faPatchList& patches = *this;

    lduInterfacePtrsList list(patches.size());

    forAll(list, patchi)
    {
        const lduInterface* lduPtr = isA<lduInterface>(patches[patchi]);

        if (lduPtr)
        {
            list.set(patchi, lduPtr);
        }
    }

    return list;
}


Foam::label Foam::faBoundaryMesh::nNonProcessor() const
{
    const faPatchList& patches = *this;

    label nonProc = 0;

    for (const faPatch& p : patches)
    {
        if (isA<processorFaPatch>(p))
        {
            break;
        }

        ++nonProc;
    }

    return nonProc;
}


const Foam::HashTable<Foam::labelList>&
Foam::faBoundaryMesh::groupPatchIDs() const
{
    if (!groupIDsPtr_)
    {
        calcGroupIDs();
    }

    return *groupIDsPtr_;
}


void Foam::faBoundaryMesh::setGroup
(
    const word& groupName,
    const labelUList& patchIDs
)
{
    groupIDsPtr_.clear();

    faPatchList& patches = *this;

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


Foam::wordList Foam::faBoundaryMesh::names() const
{
    return PtrListOps::get<word>(*this, nameOp<faPatch>());
}


Foam::wordList Foam::faBoundaryMesh::types() const
{
    return PtrListOps::get<word>(*this, typeOp<faPatch>());
}


Foam::labelList Foam::faBoundaryMesh::patchStarts() const
{
    // Manually: faPatch does not have independent start() information

    const faPatchList& patches = *this;

    labelList list(patches.size());

    label beg = mesh_.nInternalEdges();
    forAll(patches, patchi)
    {
        const label len = patches[patchi].nEdges();
        list[patchi] = beg;
        beg += len;
    }
    return list;
}


Foam::labelList Foam::faBoundaryMesh::patchSizes() const
{
    return
        PtrListOps::get<label>
        (
            *this,
            [](const faPatch& p) { return p.nEdges(); }  // avoid virtual
        );
}


Foam::List<Foam::labelRange> Foam::faBoundaryMesh::patchRanges() const
{
    const faPatchList& patches = *this;

    List<labelRange> list(patches.size());

    label beg = mesh_.nInternalEdges();
    forAll(patches, patchi)
    {
        const label len = patches[patchi].nEdges();
        list[patchi].reset(beg, len);
        beg += len;
    }
    return list;
}


Foam::label Foam::faBoundaryMesh::start() const
{
    return mesh_.nInternalEdges();
}


Foam::label Foam::faBoundaryMesh::nEdges() const
{
    return mesh_.nBoundaryEdges();
}


Foam::labelRange Foam::faBoundaryMesh::range() const
{
    return labelRange(mesh_.nInternalEdges(), mesh_.nBoundaryEdges());
}


Foam::labelList Foam::faBoundaryMesh::indices
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
        // Special version of above for reduced memory footprint

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


Foam::labelList Foam::faBoundaryMesh::indices
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


Foam::label Foam::faBoundaryMesh::findIndex(const wordRe& key) const
{
    if (key.empty())
    {
        return -1;
    }
    return PtrListOps::firstMatching(*this, key);
}


Foam::label Foam::faBoundaryMesh::findPatchID
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
            << "Available patch names: " << names() << endl
            << exit(FatalError);
    }

    // Patch not found
    if (debug)
    {
        Pout<< "label faBoundaryMesh::findPatchID(const word&) const"
            << "Patch named " << patchName << " not found.  "
            << "Available patch names: " << names() << endl;
    }

    // Not found, return -1
    return -1;
}


Foam::label Foam::faBoundaryMesh::whichPatch(const label edgeIndex) const
{
    // Find out which patch the current face belongs to by comparing label
    // with patch start labels.
    // If the face is internal, return -1;
    // if it is off the end of the list, abort
    if (edgeIndex < mesh().nInternalEdges())
    {
        return -1;
    }
    else if (edgeIndex >= mesh().nEdges())
    {
        FatalErrorInFunction
            << "Edge " << edgeIndex
            << " out of bounds. Number of geometric edges " << mesh().nEdges()
            << abort(FatalError);
    }

    forAll(*this, patchi)
    {
        label start = mesh_.patchStarts()[patchi];
        label size = operator[](patchi).faPatch::size();

        if
        (
            edgeIndex >= start
         && edgeIndex < start + size
        )
        {
            return patchi;
        }
    }

    // If not in any of above, it's trouble!
    FatalErrorInFunction
        << "error in patch search algorithm"
        << abort(FatalError);

    return -1;
}


bool Foam::faBoundaryMesh::checkParallelSync(const bool report) const
{
    if (!Pstream::parRun())
    {
        return false;
    }

    const faBoundaryMesh& bm = *this;

    bool hasError = false;

    // Collect non-proc patches and check proc patches are last.
    wordList localNames(bm.size());
    wordList localTypes(bm.size());

    label nonProci = 0;

    forAll(bm, patchi)
    {
        if (!isA<processorFaPatch>(bm[patchi]))
        {
            if (nonProci != patchi)
            {
                // A processor patch in between normal patches!
                hasError = true;

                if (debug || report)
                {
                    Pout<< " ***Problem with boundary patch " << patchi
                        << " name:" << bm[patchi].name()
                        << " type:" <<  bm[patchi].type()
                        << " - seems to be preceeded by processor patches."
                        << " This is usually a problem." << endl;
                }
            }
            else
            {
                localNames[nonProci] = bm[patchi].name();
                localTypes[nonProci] = bm[patchi].type();
                ++nonProci;
            }
        }
    }
    localNames.resize(nonProci);
    localTypes.resize(nonProci);

    // Check and report error(s) on master

    const globalIndex procAddr
    (
        // Don't need to collect master itself
        (Pstream::master() ? 0 : nonProci),
        globalIndex::gatherOnly{}
    );

    const wordList allNames(procAddr.gather(localNames));
    const wordList allTypes(procAddr.gather(localTypes));

    // Automatically restricted to master
    for (const int proci : procAddr.subProcs())
    {
        const auto procNames(allNames.slice(procAddr.range(proci)));
        const auto procTypes(allTypes.slice(procAddr.range(proci)));

        if (procNames != localNames || procTypes != localTypes)
        {
            hasError = true;

            if (debug || report)
            {
                Info<< " ***Inconsistent patches across processors, "
                       "processor0 has patch names:" << localNames
                    << " patch types:" << localTypes
                    << " processor" << proci
                    << " has patch names:" << procNames
                    << " patch types:" << procTypes
                    << endl;
            }
        }
    }

    // Reduce (not broadcast) to respect local out-of-order errors (first loop)
    reduce(hasError, orOp<bool>());

    return hasError;
}


bool Foam::faBoundaryMesh::checkDefinition(const bool report) const
{
    label nextPatchStart = mesh().nInternalEdges();
    const faBoundaryMesh& bm = *this;

    bool hasError = false;

    forAll(bm, patchi)
    {
        if (bm[patchi].start() != nextPatchStart && !hasError)
        {
            hasError = true;

            InfoInFunction
                << " ****Problem with boundary patch " << patchi
                << " named " << bm[patchi].name()
                << " of type " <<  bm[patchi].type()
                << ". The patch should start on face no " << nextPatchStart
                << " and the patch specifies " << bm[patchi].start()
                << "." << endl
                << "Possibly consecutive patches have this same problem."
                << " Suppressing future warnings." << endl;
        }

        // Warn about duplicate boundary patches?

        nextPatchStart += bm[patchi].faPatch::size();
    }

    if (hasError)
    {
        SeriousErrorInFunction
            << "This mesh is not valid: boundary definition is in error."
            << endl;
    }
    else
    {
        if (debug || report)
        {
            Info << "Boundary definition OK." << endl;
        }
    }

    return hasError;
}


void Foam::faBoundaryMesh::movePoints(const pointField& p)
{
    // processorFaPatch geometry triggers calculation of pointNormals.
    // This uses parallel comms and hence will not be trigggered
    // on processors that do not have a processorFaPatch so instead
    // force construction.
    (void)mesh_.pointAreaNormals();

    PstreamBuffers pBufs(Pstream::defaultCommsType);

    if
    (
        pBufs.commsType() == Pstream::commsTypes::blocking
     || pBufs.commsType() == Pstream::commsTypes::nonBlocking
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
    else if (pBufs.commsType() == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = mesh().globalData().patchSchedule();

        // Dummy.
        pBufs.finishedSends();

        for (const auto& schedEval : patchSchedule)
        {
            const label patchi = schedEval.patch;

            if (schedEval.init)
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


void Foam::faBoundaryMesh::updateMesh()
{
    PstreamBuffers pBufs(Pstream::defaultCommsType);

    if
    (
        pBufs.commsType() == Pstream::commsTypes::blocking
     || pBufs.commsType() == Pstream::commsTypes::nonBlocking
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
    else if (pBufs.commsType() == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = mesh().globalData().patchSchedule();

        // Dummy.
        pBufs.finishedSends();

        for (const auto& schedEval : patchSchedule)
        {
            const label patchi = schedEval.patch;

            if (schedEval.init)
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


bool Foam::faBoundaryMesh::writeData(Ostream& os) const
{
    const faPatchList& patches = *this;

    os  << patches.size() << nl << token::BEGIN_LIST << incrIndent << nl;

    for (const faPatch& p : patches)
    {
        os.beginBlock(p.name());
        os << p;
        os.endBlock();
    }

    os  << decrIndent << token::END_LIST;

    os.check(FUNCTION_NAME);
    return os.good();
}


bool Foam::faBoundaryMesh::writeObject
(
    IOstreamOption streamOpt,
    const bool valid
) const
{
    // Allow/disallow compression?
    // 1. keep readable
    // 2. save some space
    // ??? streamOpt.compression(IOstreamOption::UNCOMPRESSED);
    return regIOobject::writeObject(streamOpt, valid);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const faBoundaryMesh& bm)
{
    bm.writeData(os);
    return os;
}


// ************************************************************************* //
