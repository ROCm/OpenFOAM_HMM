/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

#include "faBoundaryMesh.H"
#include "faMesh.H"
#include "primitiveMesh.H"
#include "PtrListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faBoundaryMesh, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// bool Foam::faBoundaryMesh::hasGroupIDs() const
// {
//     /// if (groupIDsPtr_)
//     /// {
//     ///     // Use existing cache
//     ///     return !groupIDsPtr_->empty();
//     /// }
//
//     const faPatchList& patches = *this;
//
//     for (const faPatch& p : patches)
//     {
//         if (!p.inGroups().empty())
//         {
//             return true;
//         }
//     }
//
//     return false;
// }


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
    if
    (
        readOpt() == IOobject::MUST_READ
     || readOpt() == IOobject::MUST_READ_IF_MODIFIED
    )
    {
        // Warn for MUST_READ_IF_MODIFIED
        warnNoRereading<faBoundaryMesh>();

        faPatchList& patches = *this;

        // Read faPatch list
        Istream& is = readStream(typeName);

        PtrList<entry> patchEntries(is);
        patches.setSize(patchEntries.size());

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
    }
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

    forAll(*this, patchi)
    {
        operator[](patchi).initGeometry();
    }

    forAll(*this, patchi)
    {
        operator[](patchi).calcGeometry();
    }
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
    const bool useGroups  /* ignored */
) const
{
    if (matcher.empty())
    {
        return labelList();
    }

    if (matcher.isPattern())
    {
        return PtrListOps::findMatching(*this, matcher);
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
    }

    return labelList();
}


Foam::labelList Foam::faBoundaryMesh::indices
(
    const wordRes& matcher,
    const bool useGroups  /* ignored */
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

    return PtrListOps::findMatching(*this, matcher);
}


Foam::label Foam::faBoundaryMesh::findIndex(const wordRe& key) const
{
    if (key.empty())
    {
        return -1;
    }
    return PtrListOps::firstMatching(*this, key);
}


Foam::label Foam::faBoundaryMesh::findPatchID(const word& patchName) const
{
    if (patchName.empty())
    {
        return -1;
    }

    return PtrListOps::firstMatching(*this, patchName);
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

    faPatchList& patches = *this;

    forAll(patches, patchi)
    {
        patches[patchi].initMovePoints(p);
    }

    forAll(patches, patchi)
    {
        patches[patchi].movePoints(p);
    }
}


void Foam::faBoundaryMesh::updateMesh()
{
    faPatchList& patches = *this;

    forAll(patches, patchi)
    {
        patches[patchi].initUpdateMesh();
    }

    forAll(patches, patchi)
    {
        patches[patchi].updateMesh();
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
