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

#include "topoBitSet.H"
#include "polyMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::topoBitSet::topoBitSet
(
    const polyMesh& mesh,
    const word& setName
)
:
    topoSet
    (
        IOobject
        (
            setName,
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        0  // zero-sized (unallocated) labelHashSet
    ),
    selected_()
{}


Foam::topoBitSet::topoBitSet
(
    const polyMesh& mesh,
    const word& setName,
    const label size,
    const bool val
)
:
    topoBitSet(mesh, setName)
{
    selected_.resize(size, val);
}


Foam::topoBitSet::topoBitSet
(
    const polyMesh& mesh,
    const word& setName,
    const label size,
    const bitSet& bits
)
:
    topoBitSet(mesh, setName)
{
    selected_ = bits;
    selected_.resize(size);
}


Foam::topoBitSet::topoBitSet
(
    const polyMesh& mesh,
    const word& setName,
    const label size,
    bitSet&& bits
)
:
    topoBitSet(mesh, setName)
{
    selected_ = std::move(bits);
    selected_.resize(size);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::topoBitSet::found(const label id) const
{
    return selected_.test(id);
}


bool Foam::topoBitSet::set(const label id)
{
    return selected_.set(id);
}


bool Foam::topoBitSet::unset(const label id)
{
    return selected_.unset(id);
}


void Foam::topoBitSet::set(const labelUList& labels)
{
    selected_.set(labels);
}


void Foam::topoBitSet::unset(const labelUList& labels)
{
    selected_.unset(labels);
}


void Foam::topoBitSet::invert(const label maxLen)
{
    selected_.resize(maxLen);
    selected_.flip();
}


void Foam::topoBitSet::subset(const topoSet& set)
{
    // Only retain entries found in both sets
    if (isA<topoBitSet>(set))
    {
        selected_ &= refCast<const topoBitSet>(set).selected_;
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


void Foam::topoBitSet::addSet(const topoSet& set)
{
    // Add entries to the set
    if (isA<topoBitSet>(set))
    {
        selected_ |= refCast<const topoBitSet>(set).selected_;
    }
    else
    {
        for (const label id : set)
        {
            selected_.set(id);
        }
    }
}


void Foam::topoBitSet::subtractSet(const topoSet& set)
{
    // Subtract entries from the set
    if (isA<topoBitSet>(set))
    {
        selected_ -= refCast<const topoBitSet>(set).selected_;
    }
    else
    {
        for (const label id : set)
        {
            selected_.unset(id);
        }
    }
}


// ************************************************************************* //
