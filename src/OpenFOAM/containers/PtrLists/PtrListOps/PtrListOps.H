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

InNamespace
    Foam

Description
    Functions to operate on Pointer Lists.

Namespace
    Foam::PtrListOps

Description
    Various utility functions to operate on Pointer Lists.

SourceFiles
    PtrListOpsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef PtrListOps_H
#define PtrListOps_H

#include "PtrList.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Generate (stable) sort order for the list
template<class T>
void sortedOrder(const UPtrList<T>& input, labelList& order);

//- Generate (stable) sort order for the list,
//- using the specified list compare predicate
template<class T, class ListComparePredicate>
void sortedOrder
(
    const UPtrList<T>& input,
    labelList& order,
    const ListComparePredicate& comp
);


//- Inplace (stable) sorting of pointer list.
template<class T>
void sort(UPtrList<T>& list);

//- Inplace (stable) sorting of pointer list.
template<class T, class Compare>
void sort(UPtrList<T>& list, const Compare& comp);

//- Inplace shuffle of pointer list.
template<class T>
void shuffle(UPtrList<T>& list);


/*---------------------------------------------------------------------------*\
                       Namespace PtrListOps Declaration
\*---------------------------------------------------------------------------*/

namespace PtrListOps
{

// Public Classes

//- A UPtrList compare binary predicate for normal sort.
//  Null entries sort to the end
template<class T>
struct less
{
    const UPtrList<T>& values;

    less(const UPtrList<T>& list)
    :
        values(list)
    {}

    bool operator()(const label a, const label b) const
    {
        const T* const ptr1 = values(a);
        const T* const ptr2 = values(b);

        return (ptr1 && ptr2) ? (*ptr1 < *ptr2) : !ptr2;
    }
};


//- A UPtrList compare binary predicate for reverse sort.
//  Null entries sort to the end
template<class T>
struct greater
{
    const UPtrList<T>& values;

    greater(const UPtrList<T>& list)
    :
        values(list)
    {}

    bool operator()(const label a, const label b) const
    {
        const T* const ptr1 = values(a);
        const T* const ptr2 = values(b);

        return (ptr1 && ptr2) ? (*ptr2 < *ptr1) : !ptr2;
    }
};


} // End namespace ListOps

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PtrListOpsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
