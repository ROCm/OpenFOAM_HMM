/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "fvMesh.H"
#include "fvBoundaryMesh.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvBoundaryMesh::addPatches(const polyBoundaryMesh& basicBdry)
{
    setSize(basicBdry.size());

    // Set boundary patches
    fvPatchList& Patches = *this;

    forAll(Patches, patchI)
    {
        Patches.set(patchI, fvPatch::New(basicBdry[patchI], *this));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvBoundaryMesh::fvBoundaryMesh
(
    const fvMesh& m
)
:
    fvPatchList(0),
    mesh_(m)
{}


Foam::fvBoundaryMesh::fvBoundaryMesh
(
    const fvMesh& m,
    const polyBoundaryMesh& basicBdry
)
:
    fvPatchList(basicBdry.size()),
    mesh_(m)
{
    addPatches(basicBdry);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvBoundaryMesh::movePoints()
{
    forAll(*this, patchi)
    {
        operator[](patchi).initMovePoints();
    }

    forAll(*this, patchi)
    {
        operator[](patchi).movePoints();
    }
}


Foam::lduInterfacePtrsList Foam::fvBoundaryMesh::interfaces() const
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


void Foam::fvBoundaryMesh::readUpdate(const polyBoundaryMesh& basicBdry)
{
    clear();
    addPatches(basicBdry);
}


// ************************************************************************* //
