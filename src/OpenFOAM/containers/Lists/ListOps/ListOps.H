/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

InNamspace
    Foam

Description
    Various functions to operate on Lists.

SourceFiles
    ListOps.C
    ListOpsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef ListOps_H
#define ListOps_H

#include "FlatOutput.H"
#include "labelList.H"
#include "ops.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

extern const labelList emptyLabelList;

//- Return reference to zero-sized list. Compare to List::null() which returns
//  null pointer cast as list reference.
template<class Type>
static const List<Type>& emptyList()
{
    return *reinterpret_cast<const List<Type>*>(&emptyLabelList);
}


//- Renumber the values (not the indices) of a list.
//  Negative ListType elements are left untouched.
template<class ListType>
ListType renumber(const labelUList& oldToNew, const ListType& lst);

//- Inplace renumber the values (not the indices) of a list.
//  Negative ListType elements are left untouched.
template<class ListType>
void inplaceRenumber(const labelUList& oldToNew, ListType& lst);


//- Reorder the elements of a list.
//  Locations with negative oldToNew values are left as is (copy-through).
//  However, if pruning is activated, these negative oldToNew values are
//  instead skipped over and the resulting list shrunk to the max index
//  actually used.
template<class ListType>
ListType reorder
(
    const labelUList& oldToNew,
    const ListType& lst,
    const bool prune = false
);

//- Inplace reorder the elements of a list.
//  Locations with negative oldToNew values are left as is (copy-through).
//  However, if pruning is activated, these negative oldToNew values are
//  instead skipped over and the resulting list shrunk to the max index
//  actually used.
template<class ListType>
void inplaceReorder
(
    const labelUList& oldToNew,
    ListType& lst,
    const bool prune = false
);


// Variants to work with iterators and sparse tables.
// Need to have iterators and insert()

//- Map values. Do not map negative values.
template<class Container>
void inplaceMapValue(const labelUList& oldToNew, Container& lst);

//- Recreate with mapped keys. Do not map elements with negative key.
template<class Container>
void inplaceMapKey(const labelUList& oldToNew, Container& lst);


//- Generate the (stable) sort order for the list
template<class T>
void sortedOrder(const UList<T>& lst, labelList& order);

//- Sort using specified list compare predicate
template<class T, class ListComparePredicate>
void sortedOrder
(
    const UList<T>& lst,
    labelList& order,
    const ListComparePredicate& comp
);


//- Generate (sorted) indices corresponding to duplicate list values
template<class T>
void duplicateOrder(const UList<T>& lst, labelList& order);

//- Generate (sorted) indices corresponding to duplicate list values
//  sort using specified list compare predicate
template<class T, class ListComparePredicate>
void duplicateOrder
(
    const UList<T>& lst,
    labelList& order,
    const ListComparePredicate& comp
);


//- Generate (sorted) indices corresponding to unique list values
template<class T>
void uniqueOrder(const UList<T>& lst, labelList& order);

//- Generate (sorted) indices corresponding to unique list values
//  sort using specified list compare predicate
template<class T, class ListComparePredicate>
void uniqueOrder
(
    const UList<T>& lst,
    labelList& order,
    const ListComparePredicate& comp
);


//- Inplace sorting and removal of duplicates.
//  Do not use FixedList for the input list, since it doesn't resize.
template<class ListType>
void inplaceUniqueSort(ListType& lst);

//- Inplace sorting and removal of duplicates.
//  Do not use FixedList for the input list, since it doesn't resize.
template<class ListType, class ListComparePredicate>
void inplaceUniqueSort
(
    ListType& lst,
    const ListComparePredicate& comp
);


//- Extract elements of List when select is a certain value.
//  \deprecated use subsetList instead (deprecated Mar 2017)
template<class T, class ListType>
ListType subset(const UList<T>& select, const T& value, const ListType&);

//- Inplace extract elements of List when select is a certain value.
//  \deprecated use inplaceSubsetList instead (deprecated Mar 2017)
template<class T, class ListType>
void inplaceSubset(const UList<T>& select, const T& value, ListType&);


//- Extract elements of List when select is true
//  eg, to extract all selected elements:
//    subset<boolList, labelList>(selectedElems, lst);
//  Note a labelHashSet can also be used as the bool-list.
//  Do not use FixedList for the input list, since it doesn't resize.
template<class BoolListType, class ListType>
ListType subset(const BoolListType& select, const ListType& lst);

//- Inplace extract elements of List when select is true
//  eg, to extract all selected elements:
//    inplaceSubset<boolList, labelList>(selectedElems, lst);
//  Note a labelHashSet can also be used as the bool-list.
//  Do not use FixedList for the input list, since it doesn't resize.
template<class BoolListType, class ListType>
void inplaceSubset(const BoolListType& select, ListType& lst);


//- Copy a subset of the input list when predicate is true.
//  Do not use FixedList for the input list, since it doesn't resize.
template<class ListType, class UnaryPredicate>
ListType subsetList
(
    const ListType& input,
    const UnaryPredicate& pred
);

//- Inplace subset of the list when predicate is true.
//  Do not use FixedList for the input list, since it doesn't resize.
template<class ListType, class UnaryPredicate>
void inplaceSubsetList
(
    ListType& input,
    const UnaryPredicate& pred
);


//- Invert one-to-one map. Unmapped elements will be -1.
labelList invert(const label len, const labelUList& map);

//- Invert one-to-many map. Unmapped elements will be size 0.
labelListList invertOneToMany(const label len, const labelUList& map);

//- Invert many-to-many.
//  Input and output types must be inherited from List and also
//  contain ints/labels. Used, for example, for faces to pointFaces.
template<class InputIntListType, class OutputIntListType>
void invertManyToMany
(
    const label len,
    const UList<InputIntListType>& input,
    List<OutputIntListType>& output
);

template<class InputIntListType, class OutputIntListType>
List<OutputIntListType> invertManyToMany
(
    const label len,
    const UList<InputIntListType>& input
)
{
    List<OutputIntListType> output;
    invertManyToMany<InputIntListType,OutputIntListType>(len, input, output);
    return output;
}

//- Create identity map of the given length with (map[i] == i)
labelList identity(const label len);

//- Find first occurence of given element and return index,
//  return -1 if not found. Linear search.
//  \deprecated prefer UList find/found methods (deprecated Oct 2017)
template<class ListType>
label findIndex
(
    const ListType&,
    typename ListType::const_reference,
    const label start=0
);

//- Find all occurences of given element. Linear search.
template<class ListType>
labelList findIndices
(
    const ListType& l,
    typename ListType::const_reference t,
    const label start=0
);

//- Opposite of findIndices: set values at indices to given value
template<class ListType>
void setValues
(
    ListType& l,
    const labelUList& indices,
    typename ListType::const_reference t
);

//- Opposite of findIndices: set values at indices to given value
template<class ListType>
ListType createWithValues
(
    const label sz,
    typename ListType::const_reference initValue,
    const labelUList& indices,
    typename ListType::const_reference setValue
);

//- Find index of max element (and larger than given element).
//  return -1 if not found. Linear search.
template<class ListType>
label findMax(const ListType& l, const label start=0);


//- Find index of min element (and less than given element).
//  return -1 if not found. Linear search.
template<class ListType>
label findMin(const ListType& l, const label start=0);


//- Find first occurrence of given element in sorted list and return index,
//  return -1 if not found. Binary search.
template<class ListType>
label findSortedIndex
(
    const ListType& l,
    typename ListType::const_reference t,
    const label start=0
);


//- Find last element < given value in sorted list and return index,
//  return -1 if not found. Binary search.
template<class ListType, class BinaryOp>
label findLower
(
    const ListType& l,
    typename ListType::const_reference t,
    const label start,
    const BinaryOp& bop
);


//- Find last element < given value in sorted list and return index,
//  return -1 if not found. Binary search.
template<class ListType>
label findLower
(
    const ListType& l,
    typename ListType::const_reference t,
    const label start=0
);


//- To construct a List from a C array. Has extra Container type
//  to initialise e.g. wordList from arrays of char*.
//
//  \deprecated can often use initializer_list instead (deprecated Mar 2017)
template<class Container, class T, int mRows>
List<Container> initList(const T[mRows]);


//- To construct a (square) ListList from a C array. Has extra Container type
//  to initialise e.g. faceList from arrays of labels.
//
//  \deprecated can often use initializer_list instead (deprecated Mar 2017)
template<class Container, class T, int mRows, int nColumns>
List<Container> initListList(const T[mRows][nColumns]);


//- Helper class for list to append y onto the end of x
template<class T>
class ListAppendEqOp
{
public:
    void operator()(List<T>& x, const List<T>& y) const;
};


//- Helper class for list to append unique elements of y onto the end of x
template<class T>
class ListUniqueEqOp
{
public:
    void operator()(List<T>& x, const List<T>& y) const;
};


//- Reverse a list. First element becomes last element etc.
template<class ListType>
ListType reverseList(const ListType& list);


//- Inplace reversal of a list using Swap.
template<class ListType>
void inplaceReverseList(ListType& list);


//- Rotate a list by n places. If n is positive rotate clockwise/right/down.
//  If n is negative rotate anti-clockwise/left/up.
template<class ListType>
ListType rotateList(const ListType& list, const label n);


//- Inplace reversal of a list using the Reversal Block Swapping algorithm.
template<template<typename> class ListType, class DataType>
void inplaceRotateList(ListType<DataType>& list, label n);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ListOpsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
