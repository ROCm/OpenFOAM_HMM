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

#include "HashOps.H"
#include "PackedBoolList.H"

#include <algorithm>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::labelHashSet Foam::HashSetOps::used(const PackedBoolList& select)
{
    const label count = select.count();
    const label len = select.size();

    labelHashSet output(2*count);

    label used = 0;
    for (label i = 0; i < len && used < count; ++i)
    {
        if (select[i])
        {
            output.insert(i);
            ++used;
        }
    }

    return output;
}


Foam::labelHashSet Foam::HashSetOps::used(const UList<bool>& select)
{
    // We have no estimate of the size/sparsity, just assume 1/10

    const label len = select.size();

    labelHashSet output(len/10);

    for (label i = 0; i < len; ++i)
    {
        if (select[i])
        {
            output.insert(i);
        }
    }

    return output;
}


Foam::PackedBoolList Foam::HashSetOps::bitset(const labelHashSet& labels)
{
    auto const max = std::max_element(labels.cbegin(), labels.cend());
    const label len = (max.found() ? (1 + *max) : 0);

    if (len <= 0)
    {
        return PackedBoolList();
    }

    PackedBoolList output(len);

    for (const label i : labels)
    {
        if (i >= 0)
        {
            output.set(i);
        }
    }

    return output;
}


Foam::List<bool> Foam::HashSetOps::bools(const labelHashSet& labels)
{
    auto const max = std::max_element(labels.cbegin(), labels.cend());
    const label len = (max.found() ? (1 + *max) : 0);

    if (len <= 0)
    {
        return List<bool>();
    }

    List<bool> output(len, false);

    for (const label i : labels)
    {
        if (i >= 0)
        {
            output[i] = true;
        }
    }

    return output;
}


// ************************************************************************* //
