/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::SortableList

Description
    A list that is sorted upon construction or when explicitly requested
    with the sort() method.

    Uses the Foam::stableSort() algorithm.

SourceFiles
    SortableList.C

\*---------------------------------------------------------------------------*/

#ifndef SortableList_H
#define SortableList_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class SortableList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class SortableList
:
    public List<T>
{
    // Private data

        //- Original indices
        labelList indices_;


public:

    // Constructors

        //- Null constructor, sort later (eg, after assignment or transfer)
        inline SortableList();

        //- Construct given size, sort later.
        //  The indices remain empty until the list is sorted
        explicit inline SortableList(const label size);

        //- Construct given size and initial value. Sort later on
        //  The indices remain empty until the list is sorted
        inline SortableList(const label size, const T& val);

        //- Copy construct
        inline SortableList(const SortableList<T>& lst);

        //- Move construct
        inline SortableList(SortableList<T>&& lst);

        //- Copy construct from UList, sorting immediately
        explicit SortableList(const UList<T>& values);

        //- Move construct from List, sorting immediately
        SortableList(List<T>&& values);

        //- Construct given begin/end iterators.
        //  Uses std::distance to determine the size.
        template<class InputIterator>
        inline SortableList(InputIterator begIter, InputIterator endIter);

        //- Construct from an initializer list, sorting immediately
        SortableList(std::initializer_list<T> values);

        //- Construct from transferred List, sorting immediately
        explicit SortableList(const Xfer<List<T>>& values);


    // Member Functions

        //- Return the list of sorted indices. Updated every sort
        const labelList& indices() const
        {
            return indices_;
        }

        //- Return non-const access to the sorted indices. Updated every sort
        labelList& indices()
        {
            return indices_;
        }

        //- Clear the list and the indices
        void clear();

        //- Clear the indices and return a reference to the underlying List
        List<T>& shrink();

        //- (stable) sort the list (if changed after construction time).
        //  Resizes the indices as required
        void sort();

        //- Reverse (stable) sort the list
        void reverseSort();

        //- Swap content with another SortableList in constant time
        inline void swap(SortableList<T>& lst);

        //- Transfer contents to the Xfer container as a plain List
        inline Xfer<List<T>> xfer();


    // Member Operators

        //- Assignment of all entries to the given value, removing indices.
        inline void operator=(const T& val);

        //- Assignment to UList operator, removing indices. Takes linear time
        inline void operator=(const UList<T>& lst);

        //- Assignment operator. Takes linear time
        inline void operator=(const SortableList<T>& lst);

        //- Move assignment, removing indices. Takes linear time
        inline void operator=(List<T>&& lst);

        //- Move operator. Takes linear time
        inline void operator=(SortableList<T>&& lst);

        //- Assignment to an initializer list
        void operator=(std::initializer_list<T> lst);
};


// Global Functions

// Exchange contents of lists - see SortableList::swap().
template<class T>
inline void Swap(SortableList<T>& a, SortableList<T>& b);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SortableList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
