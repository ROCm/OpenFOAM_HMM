/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
Foam::labelList Foam::sortedOrder
(
    const UPtrList<T>& input
)
{
    labelList order;
    sortedOrder(input, order, typename PtrListOps::less<T>(input));
    return order;
}


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
    // List lengths must be identical. Old content is overwritten
    order.resize_nocopy(input.size());

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
    labelList order(identity(list.size()));
    Foam::shuffle(order);
    list.sortOrder(order, false);  // false = allow nullptr
}



// Templated implementation for types(), names(), etc - file-scope
template<class ReturnType, class T, class AccessOp>
Foam::List<ReturnType> Foam::PtrListOps::get
(
    const UPtrList<T>& list,
    const AccessOp& aop
)
{
    const label len = list.size();

    List<ReturnType> output(len);

    label count = 0;
    for (label i = 0; i < len; ++i)
    {
        const T* ptr = list.get(i);

        if (bool(ptr))
        {
            output[count++] = aop(*ptr);
        }
    }

    output.resize(count);

    return output;
}


template<class T, class UnaryMatchPredicate>
Foam::List<Foam::word> Foam::PtrListOps::names
(
    const UPtrList<T>& list,
    const UnaryMatchPredicate& matcher
)
{
    // Possible: const auto aop = nameOp<T>();

    const label len = list.size();

    List<word> output(len);

    label count = 0;
    for (label i = 0; i < len; ++i)
    {
        const T* ptr = list.get(i);

        if (bool(ptr))
        {
            if (matcher(ptr->name()))
            {
                output[count++] = (ptr->name());
            }
        }
    }

    output.resize(count);

    return output;
}


template<class T>
Foam::List<Foam::word> Foam::PtrListOps::names
(
    const UPtrList<T>& list
)
{
    return PtrListOps::names(list, predicates::always());
}


template<class T, class UnaryMatchPredicate>
Foam::label Foam::PtrListOps::firstMatching
(
    const UPtrList<T>& list,
    const UnaryMatchPredicate& matcher
)
{
    const label len = list.size();

    for (label i = 0; i < len; ++i)
    {
        const T* ptr = list.get(i);

        if (bool(ptr) && matcher(ptr->name()))
        {
            return i;
        }
    }

    return -1;
}


template<class T, class UnaryMatchPredicate>
Foam::labelList Foam::PtrListOps::findMatching
(
    const UPtrList<T>& list,
    const UnaryMatchPredicate& matcher
)
{
    const label len = list.size();

    labelList output(len);

    label count = 0;
    for (label i = 0; i < len; ++i)
    {
        const T* ptr = list.get(i);

        if (bool(ptr) && matcher(ptr->name()))
        {
            output[count++] = i;
        }
    }

    output.resize(count);

    return output;
}


// ************************************************************************* //
