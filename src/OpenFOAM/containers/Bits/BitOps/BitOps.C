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

#include "BitOps.H"
#include "bitSet.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
