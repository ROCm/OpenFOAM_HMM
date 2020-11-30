/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "bitSet.H"
#include "labelRange.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bitSet, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::bitSet& Foam::bitSet::minusEq(const bitSet& other)
{
    if (&other == this)
    {
        // Self '-=' : results in clearing all bits
        if (debug & 2)
        {
            InfoInFunction
                << "Perform -= on self: clears all bits" << nl;
        }

        reset();
        return *this;
    }
    else if (empty() || other.empty())
    {
        return *this;
    }


    // The operation (on overlapping blocks)
    {
        const label nblocks = num_blocks(std::min(size(), other.size()));
        const auto& rhs = other.blocks_;

        for (label blocki = 0; blocki < nblocks; ++blocki)
        {
            blocks_[blocki] &= ~rhs[blocki];
        }
    }

    return *this;
}


Foam::bitSet& Foam::bitSet::andEq(const bitSet& other)
{
    if (&other == this)
    {
        // Self '&=' : no-op

        if (debug & 2)
        {
            InfoInFunction
                << "Perform &= on self: ignore" << nl;
        }

        return *this;
    }
    else if (empty())
    {
        // empty set : no-op (no overlap possible)
        return *this;
    }
    else if (other.empty())
    {
        reset();  // Other is empty - no overlap possible
        return *this;
    }


    // The operation (on overlapping blocks)
    {
        const label nblocks = num_blocks(std::min(size(), other.size()));
        const auto& rhs = other.blocks_;

        for (label blocki = 0; blocki < nblocks; ++blocki)
        {
            blocks_[blocki] &= rhs[blocki];
        }
    }

    return *this;
}


Foam::bitSet& Foam::bitSet::orEq(const bitSet& other, const bool strict)
{
    if (&other == this)
    {
        // Self '|=' : no-op

        if (debug & 2)
        {
            InfoInFunction
                << "Perform |= on self: ignore" << nl;
        }

        return *this;
    }
    else if (other.empty())
    {
        if ((debug & 2) && !empty())
        {
            // OK if both are empty
            InfoInFunction
                << "Perform |= using empty operand: ignore" << nl;
        }

        // No (normal) overlap: no-op
        return *this;
    }
    else if (empty())
    {
        if (debug & 2)
        {
            InfoInFunction
                << "Perform |= on empty bitSet" << nl;
        }

        if (strict)
        {
            // No (normal) overlap: no-op
            return *this;
        }
    }
    else if ((debug & 2) && (size() != other.size()))
    {
        InfoInFunction
            << "Perform |= on dissimilar sized bitSets: "
            << size()  << " vs. " << other.size() << nl;
    }

    label minpos = -1; // Min trim point

    if ((size() < other.size()) && !strict)
    {
        // The size (B > A) and we are non-strict (greedy), which means we may
        // acquire additional bits from B. However, we would like to avoid
        // spurious changes in the size of A (ie, B is longer but the extra
        // bits are unset and thus don't affect the logical result).

        minpos = size();
        resize(other.size());   // Blocks now overlap
    }


    // The operation (on overlapping blocks)
    {
        const label nblocks = num_blocks(std::min(size(), other.size()));
        const auto& rhs = other.blocks_;

        for (label blocki = 0; blocki < nblocks; ++blocki)
        {
            blocks_[blocki] |= rhs[blocki];
        }
    }


    // Cleanup - minpos >= 0 means we need to check/adjust the trim point
    if (minpos >= 0)
    {
        trim(minpos); // Adjust the trim point (size)
    }
    else
    {
        clear_trailing_bits();
    }

    return *this;
}


Foam::bitSet& Foam::bitSet::xorEq(const bitSet& other, const bool strict)
{
    if (&other == this)
    {
        // Self '^=' : results in clearing all bits

        if (debug & 2)
        {
            InfoInFunction
                << "Perform ^= on self: clears all bits" << nl;
        }

        reset();
        return *this;
    }
    else if (other.empty())
    {
        if ((debug & 2) && !empty())
        {
            // OK if both are empty
            InfoInFunction
                << "Perform ^= using empty operand: ignore" << nl;
        }

        // No (normal) overlap: no-op
        return *this;
    }
    else if (empty())
    {
        if (debug & 2)
        {
            InfoInFunction
                << "Perform ^= on empty bitSet" << nl;
        }

        if (strict)
        {
            // No (normal) overlap: no-op
            return *this;
        }
    }
    else if ((debug & 2) && (size() != other.size()))
    {
        InfoInFunction
            << "Perform ^= on dissimilar sized bitSets: "
            << size()  << " vs. " << other.size() << nl;
    }

    label minpos = -1; // Min trim point

    if ((size() < other.size()) && !strict)
    {
        minpos = size();        // This logic is explained in the orEq() method
        resize(other.size());   // Blocks now overlap
    }


    // The operation (on overlapping blocks)
    {
        const label nblocks = num_blocks(std::min(size(), other.size()));
        const auto& rhs = other.blocks_;

        for (label blocki = 0; blocki < nblocks; ++blocki)
        {
            blocks_[blocki] ^= rhs[blocki];
        }
    }


    // Cleanup - minpos >= 0 means we need to check/adjust the trim point
    if (minpos >= 0)
    {
        trim(minpos); // Adjust the trim point (size)
    }
    else
    {
        clear_trailing_bits();
    }

    return *this;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bitSet::bitSet(Istream& is)
:
    PackedList<1>()
{
    is  >> *this;
}


Foam::bitSet::bitSet(const bitSet& bitset, const labelUList& addr)
:
    bitSet(addr.size())
{
    const label len = addr.size();

    for (label i = 0; i < len; ++i)
    {
        set(i, bitset.get(addr[i]));
    }
}


Foam::bitSet::bitSet(const bitSet& bitset, const labelRange& range)
:
    bitSet(range.size())
{
    label pos = range.start();
    const label len = range.size();

    for (label i = 0; i < len; ++i)
    {
        set(i, bitset.get(pos));
        ++pos;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::bitSet::assign(const UList<bool>& bools)
{
    const label len = bools.size();

    clear();
    resize(len);

    // Could also handle block-wise (in the future?)

    // Set according to indices that are true.
    for (label i = 0; i < len; ++i)
    {
        if (bools[i])
        {
            set(i);
        }
    }
}


bool Foam::bitSet::intersects(const bitSet& other) const
{
    if (size() && other.size())
    {
        const label nblocks = num_blocks(std::min(size(), other.size()));
        const auto& rhs = other.blocks_;

        for (label blocki = 0; blocki < nblocks; ++blocki)
        {
            if (bool(blocks_[blocki] & rhs[blocki]))
            {
                return true;
            }
        }
    }

    return false;
}


void Foam::bitSet::set(const labelRange& range)
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
    // - flood-fill the reset, which now corresponds to the full range.
    //
    // NB: use labelRange after() for the exclusive end-value, which
    // corresponds to our new set size.
    if (slice.after() >= size())
    {
        resize(slice.start(), false);
        resize(slice.after(), true);
        return;
    }

    // The more difficult case - everything in between.
    // 1. sequence may begin/end in the same block
    // 2. Cover more than one block
    //    a. with partial coverage in the first block
    //    b. with partial coverage in the end block

    // The begin block/offset
    unsigned int bblock = slice.first() / elem_per_block;
    unsigned int bmask  = slice.first() % elem_per_block;

    // The end block/offset
    unsigned int eblock = slice.after() / elem_per_block;
    unsigned int emask  = slice.after() % elem_per_block;

    // Transform offsets to lower bit masks
    if (bmask) bmask = mask_lower(bmask);
    if (emask) emask = mask_lower(emask);

    if (bblock == eblock)
    {
        // Same block - flll between the begin/end bits.
        // Example:
        // bmask = 0000000000001111  (lower bits)
        // emask = 0000111111111111  (lower bits)
        // -> set  0000111111110000  (xor)

        blocks_[bblock] |= (emask^bmask);
    }
    else
    {
        if (bmask)
        {
            // The first (partial) block
            // - set everything above the bmask.
            blocks_[bblock] |= (~bmask);
            ++bblock;
        }

        // Fill these blocks
        for (unsigned blocki = bblock; blocki < eblock; ++blocki)
        {
            blocks_[blocki] = (~0u);
        }

        if (emask)
        {
            // The last (partial) block.
            // - set everything below emask.
            blocks_[eblock] |= (emask);
        }
    }
}


void Foam::bitSet::unset(const labelRange& range)
{
    // Require intersection with the current bitset
    const labelRange slice = range.subset0(size());

    // Range does not intersect (invalid, empty, bitset is empty)
    if (slice.empty())
    {
        return;
    }

    // Range finishes at or beyond the right side.
    //
    // NB: use labelRange after() for the exclusive end-value, which
    // corresponds to our new set size.
    if (slice.after() >= size())
    {
        // The original size
        const label orig = size();

        resize(slice.start(), false);
        resize(orig, false);
        return;
    }


    // The more difficult case - everything in between.
    // 1. sequence may begin/end in the same block
    // 2. Cover more than one block
    //    a. with partial coverage in the first block
    //    b. with partial coverage in the end block

    // The begin block/offset
    unsigned int bblock = slice.first() / elem_per_block;
    unsigned int bmask  = slice.first() % elem_per_block;

    // The end block/offset
    unsigned int eblock = slice.after() / elem_per_block;
    unsigned int emask  = slice.after() % elem_per_block;

    // Transform offsets to lower bit masks
    if (bmask) bmask = mask_lower(bmask);
    if (emask) emask = mask_lower(emask);

    if (bblock == eblock)
    {
        // Same block - flll between the begin/end bits.
        // Example:
        // bmask = 0000000000001111  (lower bits)
        // emask = 0000111111111111  (lower bits)
        // -> set  0000111111110000  (xor)
        // -> ~    1111000000001111

        blocks_[bblock] &= (~(emask^bmask));
    }
    else
    {
        if (bmask)
        {
            // The first (partial) block
            // - only retain things below bmask.
            blocks_[bblock] &= (bmask);
            ++bblock;
        }

        // Clear these blocks
        for (unsigned blocki = bblock; blocki < eblock; ++blocki)
        {
            blocks_[blocki] = (0u);
        }

        if (emask)
        {
            // The last (partial) block.
            // - only retain things above bmask.
            blocks_[eblock] &= (~emask);
        }
    }
}


Foam::labelList Foam::bitSet::toc() const
{
    // Number of used (set) entries
    const label total = any() ? count() : 0;

    if (!total)
    {
        return labelList();
    }

    labelList output(total);
    label nItem = 0;

    // Process block-wise, detecting any '1' bits

    const label nblocks = num_blocks(size());
    for (label blocki = 0; blocki < nblocks; ++blocki)
    {
        unsigned int blockval = blocks_[blocki];

        if (blockval)
        {
            for (label pos = (blocki * elem_per_block); blockval; ++pos)
            {
                if (blockval & 1u)
                {
                    output[nItem] = pos;
                    ++nItem;
                }
                blockval >>= 1u;
            }
            if (nItem == total) break;  // Terminate early
        }
    }

    return output;
}


Foam::List<bool> Foam::bitSet::values() const
{
    List<bool> output(size(), false);

    // Process block-wise, detecting any '1' bits

    const label nblocks = num_blocks(size());
    for (label blocki = 0; blocki < nblocks; ++blocki)
    {
        label pos = (blocki * elem_per_block);

        for
        (
            unsigned int blockval = blocks_[blocki];
            blockval;
            blockval >>= 1u
        )
        {
            if (blockval & 1u)
            {
                output[pos] = true;
            }
            ++pos;
        }
    }

    return output;
}


// ************************************************************************* //
