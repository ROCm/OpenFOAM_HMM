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

#include "PackedList.H"
#include "labelRange.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<unsigned Width>
Foam::PackedList<Width>::PackedList
(
    const PackedList<Width>& list,
    const labelUList& addr
)
:
    PackedList<Width>(addr.size())
{
    const label len = addr.size();

    for (label i = 0; i < len; ++i)
    {
        set(i, list.get(addr[i]));
    }
}


template<unsigned Width>
Foam::PackedList<Width>::PackedList
(
    const PackedList<Width>& list,
    const labelUIndList& addr
)
:
    PackedList<Width>(addr.size())
{
    const label len = addr.size();

    for (label i = 0; i < len; ++i)
    {
        set(i, list.get(addr[i]));
    }
}


template<unsigned Width>
Foam::PackedList<Width>::PackedList
(
    const PackedList<Width>& list,
    const labelRange& range
)
:
    PackedList<Width>(range.size())
{
    label pos = range.start();
    const label len = range.size();

    for (label i = 0; i < len; ++i)
    {
        set(i, list.get(pos));
        ++pos;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<unsigned Width>
bool Foam::PackedList<Width>::uniform() const
{
    if (size() < 2)
    {
        return false;   // Trivial case
    }

    // The value of the first element for testing
    const unsigned int val = get(0);

    const label nblocks = num_blocks(size());

    bool identical = true;

    if (!val)
    {
        // Zero value: can just check block content directly

        for (label blocki = 0; identical && blocki < nblocks; ++blocki)
        {
            identical = !blocks_[blocki];
        }

        return identical;
    }
    else if (nblocks > 1)
    {
        // Check all blocks that are completely occupied: (nblocks-1)
        const unsigned int blockval =
            BitOps::repeat_value<block_type,Width>(val);

        for (label blocki = 0; identical && blocki < (nblocks-1); ++blocki)
        {
            identical = (blocks_[blocki] == blockval);
        }
    }

    // Partial block: check manually
    for
    (
        label elemi = elem_per_block*(nblocks-1);
        identical && elemi < size();
        ++elemi
    )
    {
        identical = (val == get(elemi));
    }

    return identical;
}


template<unsigned Width>
Foam::labelList Foam::PackedList<Width>::values() const
{
    return this->unpack<label>();
}


template<unsigned Width>
template<class IntType>
Foam::List<IntType>
Foam::PackedList<Width>::unpack() const
{
    static_assert
    (
        std::is_integral<IntType>::value,
        "Integral required for output."
    );
    static_assert
    (
        std::numeric_limits<IntType>::digits >= Width,
        "Width of IntType is too small to hold result"
    );

    if (size() < 2 || uniform())
    {
        const IntType val = (size() ? get(0) : 0);

        return List<IntType>(size(), val);
    }

    List<IntType> output(size());
    label outi = 0;

    // Process n-1 complete blocks
    const label nblocks = num_blocks(size());

    for (label blocki=0; blocki < nblocks-1; ++blocki)
    {
        unsigned int blockval = blocks_[blocki];

        for (unsigned nget = elem_per_block; nget; --nget, ++outi)
        {
            output[outi] = IntType(blockval & max_value);
            blockval >>= Width;
        }
    }

    // Any partial blocks
    for (/*nil*/; outi < size(); ++outi)
    {
        output[outi] = get(outi);
    }

    return output;
}


template<unsigned Width>
template<class IntType>
Foam::List<IntType>
Foam::PackedList<Width>::unpack(const labelRange& range) const
{
    static_assert
    (
        std::is_integral<IntType>::value,
        "Integral required for unpack output."
    );
    static_assert
    (
        std::numeric_limits<IntType>::digits >= Width,
        "Width of IntType is too small to hold unpack output."
    );


    // Could be more efficient but messier with block-wise access.
    // - automatically handles any invalid positions

    auto pos = range.start();

    List<IntType> output(range.size());

    for (IntType& out : output)
    {
        out = IntType(get(pos));
        ++pos;
    }

    return output;
}


template<unsigned Width>
template<class IntType>
Foam::List<IntType>
Foam::PackedList<Width>::unpack(const labelUList& locations) const
{
    static_assert
    (
        std::is_integral<IntType>::value,
        "Integral required for unpack output."
    );
    static_assert
    (
        std::numeric_limits<IntType>::digits >= Width,
        "Width of IntType is too small to hold unpack output."
    );


    label pos = 0;

    List<IntType> output(locations.size());

    for (IntType& out : output)
    {
        out = IntType(get(locations[pos]));
        ++pos;
    }

    return output;
}


// ************************************************************************* //
