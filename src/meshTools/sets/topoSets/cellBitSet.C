/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

#include "cellBitSet.H"
#include "polyMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellBitSet, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellBitSet::cellBitSet(const polyMesh& mesh)
:
    cellBitSet(mesh, false)
{}


Foam::cellBitSet::cellBitSet(const polyMesh& mesh, const bool val)
:
    topoSet
    (
        IOobject
        (
            "cellBitSet",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        0  // zero-sized (unallocated) labelHashSet
    ),
    selected_(mesh.nCells(), val)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cellBitSet::found(const label id) const
{
    return selected_.test(id);
}


bool Foam::cellBitSet::set(const label id)
{
    return selected_.set(id);
}


bool Foam::cellBitSet::unset(const label id)
{
    return selected_.unset(id);
}


void Foam::cellBitSet::set(const labelUList& labels)
{
    selected_.set(labels);
}


void Foam::cellBitSet::unset(const labelUList& labels)
{
    selected_.unset(labels);
}


void Foam::cellBitSet::invert(const label maxLen)
{
    selected_.resize(maxLen);
    selected_.flip();
}


void Foam::cellBitSet::subset(const topoSet& set)
{
    // Only retain entries found in both sets
    if (isA<cellBitSet>(set))
    {
        selected_ &= refCast<const cellBitSet>(set).selected_;
    }
    else if (set.empty())
    {
        selected_.reset();
    }
    else
    {
        for (const label id : selected_)
        {
            if (!set.found(id))
            {
                selected_.unset(id);
            }
        }
    }
}


void Foam::cellBitSet::addSet(const topoSet& set)
{
    // Add entries to the set
    if (isA<cellBitSet>(set))
    {
        selected_ |= refCast<const cellBitSet>(set).selected_;
    }
    else
    {
        for (const label id : set)
        {
            selected_.set(id);
        }
    }
}


void Foam::cellBitSet::subtractSet(const topoSet& set)
{
    // Subtract entries from the set
    if (isA<cellBitSet>(set))
    {
        selected_ -= refCast<const cellBitSet>(set).selected_;
    }
    else
    {
        for (const label id : set)
        {
            selected_.unset(id);
        }
    }
}


Foam::label Foam::cellBitSet::maxSize(const polyMesh& mesh) const
{
    return mesh.nCells();
}


void Foam::cellBitSet::writeDebug
(
    Ostream& os,
    const primitiveMesh& mesh,
    const label maxLen
) const
{
    topoSet::writeDebug(os, mesh.cellCentres(), maxLen);
}


// ************************************************************************* //
