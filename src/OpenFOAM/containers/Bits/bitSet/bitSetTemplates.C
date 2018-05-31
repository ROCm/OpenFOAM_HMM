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

#include <algorithm>

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class InputIter>
Foam::label Foam::bitSet::setMany(InputIter first, InputIter last)
{
    // Check the max expected value first
    const auto max = std::max_element(first, last);
    const label len = (max != last ? (1 + *max) : 0);

    label changed = 0;

    if (len > 0)
    {
        reserve(len);

        for (; first != last; ++first)
        {
            if (set(*first))
            {
                ++changed;
            }
        }
    }

    return changed;
}


template<class InputIter>
Foam::label Foam::bitSet::unset(InputIter first, InputIter last)
{
    label changed = 0;

    for (; first != last; ++first)
    {
        if (unset(*first))
        {
            ++changed;
        }
    }

    return changed;
}


// ************************************************************************* //
