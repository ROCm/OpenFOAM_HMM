/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2016-2017 Wikki Ltd
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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faBoundaryMesh, 0);
}

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    // Templated implementation for types(), names(), etc - file-scope
    template<class ListType, class GetOp>
    static ListType getMethodImpl
    (
        const faPatchList& list,
        const GetOp& getop
    )
    {
        const label len = list.size();

        ListType output(len);

        for (label i = 0; i < len; ++i)
        {
            output[i] = getop(list[i]);
        }

        return output;
    }


    // Templated implementation for indices() - file-scope
    template<class UnaryMatchPredicate>
    static labelList indicesImpl
    (
        const faPatchList& list,
        const UnaryMatchPredicate& matcher
    )
    {
        const label len = list.size();

        labelList output(len);

        label count = 0;
        for (label i = 0; i < len; ++i)
        {
            if (matcher(list[i].name()))
            {
                output[count++] = i;
            }
        }

        output.resize(count);

        return output;
    }


    // Templated implementation for findIndex() - file-scope
    template<class UnaryMatchPredicate>
    label findIndexImpl
    (
        const faPatchList& list,
        const UnaryMatchPredicate& matcher
    )
    {
        const label len = list.size();

        for (label i = 0; i < len; ++i)
        {
            if (matcher(list[i].name()))
            {
                return i;
            }
        }

        return -1;
    }

} // End namespace Foam


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
    if (readOpt() == IOobject::MUST_READ)
    {
        faPatchList& patches = *this;

        // Read polyPatchList
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


const Foam::faMesh& Foam::faBoundaryMesh::mesh() const
{
    return mesh_;
}


Foam::lduInterfacePtrsList Foam::faBoundaryMesh::interfaces() const
{
    lduInterfacePtrsList interfaces(size());

    forAll(interfaces, patchi)
    {
        if (isA<lduInterface>(this->operator[](patchi)))
        {
            interfaces.set
            (
                patchi,
                &refCast<const lduInterface>(this->operator[](patchi))
            );
        }
    }

    return interfaces;
}


Foam::wordList Foam::faBoundaryMesh::names() const
{
    return getMethodImpl<wordList>(*this, getNameOp<faPatch>());
}


Foam::wordList Foam::faBoundaryMesh::types() const
{
    return getMethodImpl<wordList>(*this, getTypeOp<faPatch>());
}


Foam::labelList Foam::faBoundaryMesh::indices
(
    const keyType& key,
    const bool useGroups  // ignored
) const
{
    if (key.empty())
    {
        return labelList();
    }

    if (key.isPattern())
    {
        const regExp matcher(key);

        return indicesImpl(*this, matcher);
    }
    else
    {
        // Literal string.
        // Special version of above for reduced memory footprint

        const word& matcher = key;

        const label patchId = findIndexImpl(*this, matcher);

        if (patchId >= 0)
        {
            return labelList(one(), patchId);
        }
    }

    return labelList();
}


Foam::label Foam::faBoundaryMesh::findIndex(const keyType& key) const
{
    if (key.empty())
    {
        return -1;
    }
    else if (key.isPattern())
    {
        // Find as regex
        regExp matcher(key);
        return findIndexImpl(*this, matcher);
    }
    else
    {
        // Find as literal string
        const word& matcher = key;
        return findIndexImpl(*this, matcher);
    }
}


Foam::label Foam::faBoundaryMesh::findPatchID(const word& patchName) const
{
    return findIndexImpl(*this, patchName);
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


Foam::Ostream& Foam::operator<<(Ostream& os, const faBoundaryMesh& bm)
{
    bm.writeData(os);
    return os;
}


// ************************************************************************* //
