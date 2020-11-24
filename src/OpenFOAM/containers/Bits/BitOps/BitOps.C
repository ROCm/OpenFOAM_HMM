/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "BitOps.H"
#include "bitSet.H"
#include "HashSet.H"
#include "List.H"
#include "labelRange.H"

// * * * * * * * * * * * * * * * * * BitOps  * * * * * * * * * * * * * * * * //

// See bitSet::set(labelRange) for original implementation
void Foam::BitOps::set(List<bool>& bools, const labelRange& range)
{
    labelRange slice(range);
    slice.adjust();  // No negative start, size adjusted accordingly

    // Range is invalid (zero-sized or entirely negative) - noop
    if (slice.empty())
    {
        return;
    }

    // Range finishes at or beyond the right side.
    // - zero fill any gaps that we might create.
    // - flood-fill the rest, which now corresponds to the full range.
    //
    // NB: use labelRange after() for the exclusive end-value, which
    // corresponds to our new set size.
    if (slice.after() >= bools.size())
    {
        label i = bools.size();

        bools.resize(slice.after(), true);

        // Backfill with false
        while (i < slice.start())
        {
            bools.unset(i);
            ++i;
        }
        return;
    }

    for (label i = slice.first(); i <= slice.last(); ++i)
    {
        bools.set(i);
    }
}


// See bitSet::set(labelRange) for original implementation
void Foam::BitOps::set(labelHashSet& hashset, const labelRange& range)
{
    labelRange slice(range);
    slice.adjust();  // No negative start, size adjusted accordingly

    for (label i = slice.first(); i <= slice.last(); ++i)
    {
        hashset.set(i);
    }
}


void Foam::BitOps::set(bitSet& bitset, const labelRange& range)
{
    bitset.set(range);
}


// See bitSet::unset(labelRange) for original implementation
void Foam::BitOps::unset(List<bool>& bools, const labelRange& range)
{
    for (label i = range.first(); i <= range.last(); ++i)
    {
        bools.unset(i);
    }
}


void Foam::BitOps::unset(labelHashSet& hashset, const labelRange& range)
{
    for (label i = range.first(); i <= range.last(); ++i)
    {
        hashset.unset(i);
    }
}


void Foam::BitOps::unset(bitSet& bitset, const labelRange& range)
{
    bitset.unset(range);
}


// * * * * * * * * * * * * * * * * BitSetOps * * * * * * * * * * * * * * * * //

Foam::bitSet Foam::BitSetOps::create
(
    const label n,
    const labelHashSet& locations,
    const bool on
)
{
    bitSet output(n, !on);

    for (const label idx : locations)
    {
        // Restrict the input size
        if (idx < n)
        {
            output.set(idx, on);
        }
    }

    return output;
}


Foam::bitSet Foam::BitSetOps::create
(
    const label n,
    const labelUList& locations,
    const bool on
)
{
    bitSet output(n, !on);

    for (const label idx : locations)
    {
        // Restrict the input size
        if (idx < n)
        {
            output.set(idx, on);
        }
    }

    return output;
}


Foam::bitSet Foam::BitSetOps::create
(
    const label n,
    const label select,
    const labelUList& values,
    const bool on
)
{
    bitSet output(n, !on);

    // Restrict the input size
    const label len = std::min(n, values.size());

    for (label idx = 0; idx < len; ++idx)
    {
        if (select == values[idx])
        {
            output.set(idx, on);
        }
    }

    return output;
}


// ************************************************************************* //
