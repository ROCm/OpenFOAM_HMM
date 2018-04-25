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
    if (size() < 2 || uniform())
    {
        const label val = (size() ? get(0) : 0);

        return labelList(size(), val);
    }

    labelList output(size());
    label outi = 0;

    // Process n-1 complete blocks
    const label nblocks = num_blocks(size());

    for (label blocki=0; blocki < nblocks-1; ++blocki)
    {
        unsigned int blockval = blocks_[blocki];

        for (unsigned nget = elem_per_block; nget; --nget, ++outi)
        {
            output[outi] = label(blockval & max_value);
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


// ************************************************************************* //
