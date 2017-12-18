/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
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

Description

\*---------------------------------------------------------------------------*/

#include "faBoundaryMesh.H"
#include "faMesh.H"
#include "primitiveMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faBoundaryMesh, 0);
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
    if (readOpt() == IOobject::MUST_READ)
    {
        faPatchList& patches = *this;

        // Read polyPatchList
        Istream& is = readStream(typeName);

        PtrList<entry> patchEntries(is);
        patches.setSize(patchEntries.size());

        forAll(patches, patchI)
        {
            patches.set
            (
                patchI,
                faPatch::New
                (
                    patchEntries[patchI].keyword(),
                    patchEntries[patchI].dict(),
                    patchI,
                    *this
                )
            );
        }

        // Check state of IOstream
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


Foam::wordList Foam::faBoundaryMesh::types() const
{
    const faPatchList& patches = *this;

    wordList t(patches.size());

    forAll(patches, patchI)
    {
        t[patchI] = patches[patchI].type();
    }

    return t;
}


Foam::wordList Foam::faBoundaryMesh::names() const
{
    const faPatchList& patches = *this;

    wordList t(patches.size());

    forAll(patches, patchI)
    {
        t[patchI] = patches[patchI].name();
    }

    return t;
}


Foam::label Foam::faBoundaryMesh::findPatchID(const word& patchName) const
{
    const faPatchList& patches = *this;

    forAll(patches, patchI)
    {
        if (patches[patchI].name() == patchName)
        {
            return patchI;
        }
    }

    // Patch not found
    return -1;
}


Foam::labelList Foam::faBoundaryMesh::findIndices
(
    const keyType& key,
    const bool usePatchGroups
) const
{
    DynamicList<label> indices;

    if (!key.empty())
    {
        if (key.isPattern())
        {
            indices = findStrings(key, this->names());
        }
        else
        {
            // Literal string. Special version of above to avoid
            // unnecessary memory allocations

            indices.setCapacity(1);
            forAll(*this, i)
            {
                if (key == operator[](i).name())
                {
                    indices.append(i);
                    break;
                }
            }
        }
    }

    return indices;
}


Foam::label Foam::faBoundaryMesh::whichPatch(const label edgeIndex) const
{
    // Find out which patch the current face belongs to by comparing label
    // with patch start labels.
    // If the face is internal, return -1;
    // if it is off the end of the list, abort
    if (edgeIndex >= mesh().nEdges())
    {
        FatalErrorInFunction
            << "given label greater than the number of edges"
            << abort(FatalError);
    }

    if (edgeIndex < mesh().nInternalEdges())
    {
        return -1;
    }

    forAll(*this, patchI)
    {
        label start = mesh_.patchStarts()[patchI];
        label size = operator[](patchI).faPatch::size();

        if
        (
            edgeIndex >= start
         && edgeIndex < start + size
        )
        {
            return patchI;
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

    bool boundaryError = false;

    forAll(bm, patchI)
    {
        if (bm[patchI].start() != nextPatchStart)
        {
            boundaryError = true;

            InfoInFunction
                << "Problem with boundary patch " << patchI
                << ".\nThe patch should start on face no " << nextPatchStart
                << " and the boundary file specifies " << bm[patchI].start()
                << "." << nl << endl;
        }

        nextPatchStart += bm[patchI].faPatch::size();
    }

    if (boundaryError)
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

    return boundaryError;
}


void Foam::faBoundaryMesh::movePoints(const pointField& p)
{
    faPatchList& patches = *this;

    forAll(patches, patchI)
    {
        patches[patchI].initMovePoints(p);
    }

    forAll(patches, patchI)
    {
        patches[patchI].movePoints(p);
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

    forAll(patches, patchi)
    {
        os  << indent << patches[patchi].name() << nl
            << indent << token::BEGIN_BLOCK << nl
            << incrIndent << patches[patchi] << decrIndent
            << indent << token::END_BLOCK << endl;
    }

    os  << decrIndent << token::END_LIST;

    // Check state of IOstream
    os.check(FUNCTION_NAME);

    return os.good();
}


Foam::Ostream& Foam::operator<<(Ostream& os, const faBoundaryMesh& bm)
{
    bm.writeData(os);
    return os;
}


// ************************************************************************* //
