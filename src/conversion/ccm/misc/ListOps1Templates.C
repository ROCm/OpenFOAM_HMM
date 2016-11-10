/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include "ListOps1.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class ListType>
ListType Foam::reorder
(
    const labelUList& oldToNew,
    const ListType& lst,
    const bool prune
)
{
    const label sz = lst.size();

    // Create copy
    ListType newLst(sz);

    // Ensure consistent addressable size (eg, DynamicList)
    newLst.setSize(sz);


    label maxIdx = 0;
    forAll(lst, elemI)
    {
        const label newIdx = oldToNew[elemI];
        if (newIdx >= 0) // could also require newIdx < sz
        {
            newLst[newIdx] = lst[elemI];
            if (prune && maxIdx < newIdx)
            {
                maxIdx = newIdx;
            }
        }
        else if (!prune)
        {
            newLst[elemI] = lst[elemI];
        }
    }

    if (prune && maxIdx < sz)
    {
        newLst.setSize(maxIdx);
    }

    return newLst;
}


template<class ListType>
void Foam::inplaceReorder
(
    const labelUList& oldToNew,
    ListType& lst,
    const bool prune
)
{
    const label sz = lst.size();

    // Create copy
    ListType newLst(sz);

    // Ensure consistent addressable size (eg, DynamicList)
    newLst.setSize(sz);

    label maxIdx = 0;
    forAll(lst, elemI)
    {
        const label newIdx = oldToNew[elemI];
        if (newIdx >= 0) // could also require newIdx < sz
        {
            newLst[newIdx] = lst[elemI];
            if (prune && maxIdx < newIdx)
            {
                maxIdx = newIdx;
            }
        }
        else if (!prune)
        {
            newLst[elemI] = lst[elemI];
        }
    }

    if (prune && maxIdx < sz)
    {
        newLst.setSize(maxIdx);
    }

    lst.transfer(newLst);
}


// ************************************************************************* //
