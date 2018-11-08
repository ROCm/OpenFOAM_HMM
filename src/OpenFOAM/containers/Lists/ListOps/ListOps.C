/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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

#include "ListOps.H"
#include <numeric>

// * * * * * * * * * * * * * * Global Data Members * * * * * * * * * * * * * //

const Foam::labelList Foam::emptyLabelList;


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::invert
(
    const label len,
    const labelUList& map
)
{
    labelList inverse(len, -1);

    forAll(map, i)
    {
        const label newIdx = map[i];

        if (newIdx >= 0)
        {
            if (inverse[newIdx] >= 0)
            {
                FatalErrorInFunction
                    << "Map is not one-to-one. At index " << i
                    << " element " << newIdx << " has already occurred before"
                    << nl << "Please use invertOneToMany instead"
                    << abort(FatalError);
            }

            inverse[newIdx] = i;
        }
    }

    return inverse;
}


Foam::labelListList Foam::invertOneToMany
(
    const label len,
    const labelUList& map
)
{
    labelList sizes(len, 0);

    for (const label newIdx : map)
    {
        if (newIdx >= 0)
        {
            sizes[newIdx]++;
        }
    }

    labelListList inverse(len);

    for (label i=0; i < len; ++i)
    {
        inverse[i].resize(sizes[i]);
        sizes[i] = 0; // reset size counter
    }

    forAll(map, i)
    {
        const label newIdx = map[i];

        if (newIdx >= 0)
        {
            inverse[newIdx][sizes[newIdx]++] = i;
        }
    }

    return inverse;
}


Foam::labelList Foam::identity(const label len, label start)
{
    labelList map(len);
    std::iota(map.begin(), map.end(), start);

    return map;
}


Foam::bitSet Foam::reorder
(
    const labelUList& oldToNew,
    const bitSet& input,
    const bool prune
)
{
    const label len = input.size();

    bitSet output;
    output.reserve(len);

    for
    (
        label pos = input.find_first();
        pos >= 0 && pos < len;
        pos = input.find_next(pos)
    )
    {
        const label newIdx = oldToNew[pos];

        if (newIdx >= 0)
        {
            output.set(newIdx);
        }
        else if (!prune)
        {
            output.set(pos);
        }
    }

    if (prune)
    {
        output.trim();
    }

    return output;
}


void Foam::inplaceReorder
(
    const labelUList& oldToNew,
    bitSet& input,
    const bool prune
)
{
    input = reorder(oldToNew, input, prune);
}


// ************************************************************************* //
