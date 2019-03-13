/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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

#include "PtrListOps.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::sortedOrder
(
    const UPtrList<T>& input,
    labelList& order
)
{
    sortedOrder(input, order, typename PtrListOps::less<T>(input));
}


template<class T, class ListComparePredicate>
void Foam::sortedOrder
(
    const UPtrList<T>& input,
    labelList& order,
    const ListComparePredicate& comp
)
{
    const label len = input.size();

    // List lengths must be identical
    if (order.size() != len)
    {
        // Avoid copying elements, they are overwritten anyhow
        order.clear();
        order.resize(len);
    }

    ListOps::identity(order);

    Foam::stableSort(order, comp);
}


template<class T>
void Foam::sort(UPtrList<T>& list)
{
    labelList order;
    sortedOrder(list, order);
    list.sortOrder(order, false);  // false = allow nullptr
}


template<class T, class Compare>
void Foam::sort(UPtrList<T>& list, const Compare& comp)
{
    labelList order;
    sortedOrder(list, order, comp);
    list.sortOrder(order, false);  // false = allow nullptr
}


template<class T>
void Foam::shuffle(UPtrList<T>& list)
{
    labelList order = identity(list.size());
    Foam::shuffle(order);
    list.sortOrder(order, false);  // false = allow nullptr
}


// ************************************************************************* //
