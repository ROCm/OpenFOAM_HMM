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

#include "topoBoolSet.H"
#include "polyMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::topoBoolSet::topoBoolSet
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


Foam::topoBoolSet::topoBoolSet
(
    const polyMesh& mesh,
    const word& setName,
    const label size,
    const bool val
)
:
    topoBoolSet(mesh, setName)
{
    selected_.resize(size, val);
}


Foam::topoBoolSet::topoBoolSet
(
    const polyMesh& mesh,
    const word& setName,
    const label size,
    const boolList& bools
)
:
    topoBoolSet(mesh, setName)
{
    selected_ = bools;
    selected_.resize(size);
}


Foam::topoBoolSet::topoBoolSet
(
    const polyMesh& mesh,
    const word& setName,
    const label size,
    boolList&& bools
)
:
    topoBoolSet(mesh, setName)
{
    selected_ = std::move(bools);
    selected_.resize(size);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::topoBoolSet::found(const label id) const
{
    return selected_.test(id);
}


bool Foam::topoBoolSet::set(const label id)
{
    return selected_.set(id);
}


bool Foam::topoBoolSet::unset(const label id)
{
    return selected_.unset(id);
}


void Foam::topoBoolSet::set(const labelUList& labels)
{
    for (const label id : labels)
    {
        selected_[id] = true;
    }
}


void Foam::topoBoolSet::unset(const labelUList& labels)
{
    for (const label id : labels)
    {
        selected_.unset(id);
    }
}


void Foam::topoBoolSet::invert(const label maxLen)
{
    selected_.resize(maxLen);
    for (bool& b : selected_)
    {
        b = !b;
    }
}


void Foam::topoBoolSet::subset(const topoSet& set)
{
    // Only retain entries found in both sets
    if (set.empty())
    {
        selected_ = false;
    }
    else
    {
        forAll(selected_, i)
        {
            selected_[i] = (selected_[i] && set.found(i));
        }
    }
}


void Foam::topoBoolSet::addSet(const topoSet& set)
{
    // Add entries to the set
    for (const label id : set)
    {
        selected_[id] = true;
    }
}


void Foam::topoBoolSet::subtractSet(const topoSet& set)
{
    // Subtract entries from the set
    for (const label id : set)
    {
        selected_.unset(id);
    }
}


// ************************************************************************* //
