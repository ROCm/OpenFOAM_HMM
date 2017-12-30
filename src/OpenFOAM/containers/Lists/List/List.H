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
    Foam::List

Description
    A 1D array of objects of type \<T\>, where the size of the vector
    is known and used for subscript bounds checking, etc.

    Storage is allocated on free-store during construction.

SourceFiles
    List.C
    ListI.H
    ListIO.C

\*---------------------------------------------------------------------------*/

#ifndef List_H
#define List_H

#include "UList.H"
#include "autoPtr.H"
#include "Xfer.H"
#include <initializer_list>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class Istream;
class Ostream;

// Forward declaration of friend functions and operators
template<class T> class List;

template<class T> Istream& operator>>(Istream& is, List<T>& L);

template<class T, unsigned Size> class FixedList;
template<class T> class PtrList;

class SLListBase;
template<class LListBase, class T> class LList;
template<class T>
using SLList = LList<SLListBase, T>;

template<class T, int SizeMin> class DynamicList;

template<class T> class SortableList;
template<class T> class IndirectList;
template<class T> class UIndirectList;
template<class T> class BiIndirectList;

typedef List<char> charList;

/*---------------------------------------------------------------------------*\
                           Class List Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class List
:
    public UList<T>
{
    // Private member functions

        //- Allocate list storage
        inline void alloc();

        //- Reallocate list storage to the given size
        inline void reAlloc(const label s);

        //- Copy list of given type
        template<class List2>
        inline void copyList(const List2& lst);

        //- Allocate storage and copy list of given type
        template<class List2>
        inline void allocCopyList(const List2& lst);

        //- Construct given begin/end iterators and number of elements
        //  Since the size is provided, the end iterator is actually ignored.
        template<class InputIterator>
        inline List
        (
            InputIterator begIter,
            InputIterator endIter,
            const label s
        );


public:

    // Static Member Functions

        //- Return a null List
        inline static const List<T>& null();


    // Constructors

        //- Null constructor
        inline List();

        //- Construct with given size
        explicit List(const label s);

        //- Construct with given size and value for all elements
        List(const label s, const T& val);

        //- Construct with given size initializing all elements to zero
        List(const label s, const zero);

        //- Copy constructor from list
        List(const List<T>& a);

        //- Copy constructor from list containing another type.
        //  This is primarily useful to convert a list of ints into floats,
        //  for example.
        template<class T2>
        explicit List(const List<T2>& a);

        //- Construct as copy or re-use as specified
        List(List<T>& a, bool reuse);

        //- Construct as subset
        List(const UList<T>& lst, const labelUList& mapAddressing);

        //- Construct given begin/end iterators.
        //  Uses std::distance to determine the size.
        template<class InputIterator>
        List(InputIterator begIter, InputIterator endIter);

        //- Construct as copy of FixedList<T, Size>
        template<unsigned Size>
        explicit List(const FixedList<T, Size>& lst);

        //- Construct as copy of PtrList<T>
        explicit List(const PtrList<T>& lst);

        //- Construct as copy of SLList<T>
        explicit List(const SLList<T>& lst);

        //- Construct as copy of UIndirectList<T>
        explicit List(const UIndirectList<T>& lst);

        //- Construct as copy of BiIndirectList<T>
        explicit List(const BiIndirectList<T>& lst);

        //- Construct from an initializer list
        List(std::initializer_list<T> lst);

        //- Transfer (move) construct
        List(const Xfer<List<T>>& lst);

        //- Move construct from List
        List(List<T>&& lst);

        //- Move construct from DynamicList
        template<int SizeMin>
        List(DynamicList<T, SizeMin>&& lst);

        //- Move construct from SortableList
        List(SortableList<T>&& lst);

        //- Construct from Istream
        List(Istream& is);

        //- Clone
        inline autoPtr<List<T>> clone() const;


    //- Destructor
    ~List();


    // Related types

        //- Declare type of subList
        typedef SubList<T> subList;


    // Member Functions

        // Edit

            //- Alias for setSize(const label)
            inline void resize(const label newSize);

            //- Alias for setSize(const label, const T&)
            inline void resize(const label newSize, const T& val);

            //- Reset size of List
            void setSize(const label newSize);

            //- Reset size of List and value for new elements
            void setSize(const label newSize, const T& val);

            //- Clear the list, i.e. set size to zero
            inline void clear();

            //- Append an element at the end of the list
            inline void append(const T& val);

            //- Move append an element at the end of the list
            inline void append(T&& val);

            //- Append a List to the end of this list
            inline void append(const UList<T>& lst);

            //- Append a UIndirectList at the end of this list
            inline void append(const UIndirectList<T>& lst);

            //- Transfer the contents of the argument List into this list
            //- and annul the argument list
            void transfer(List<T>& lst);

            //- Transfer the contents of the argument List into this list
            //- and annul the argument list
            template<int SizeMin>
            void transfer(DynamicList<T, SizeMin>& lst);

            //- Transfer the contents of the argument List into this list
            //- and annul the argument list
            void transfer(SortableList<T>& lst);

            //- Transfer contents to the Xfer container
            inline Xfer<List<T>> xfer();

            //- Return subscript-checked element of UList and resizing the list
            //- if required.
            inline T& newElmt(const label i);


        //- Disallow implicit shallowCopy
        void shallowCopy(const UList<T>&) = delete;


    // Member operators

        //- Assignment to UList operator. Takes linear time
        void operator=(const UList<T>& a);

        //- Assignment operator. Takes linear time
        void operator=(const List<T>& lst);

        //- Assignment to SLList operator. Takes linear time
        void operator=(const SLList<T>& lst);

        //- Assignment to UIndirectList operator. Takes linear time
        void operator=(const UIndirectList<T>& lst);

        //- Assignment to BiIndirectList operator. Takes linear time
        void operator=(const BiIndirectList<T>& lst);

        //- Assignment to an initializer list
        void operator=(std::initializer_list<T> lst);

        //- Assignment of all entries to the given value
        inline void operator=(const T& val);

        //- Assignment of all entries to zero
        inline void operator=(const zero);

        //- Move assignment. Takes constant time
        void operator=(List<T>&& lst);

        //- Move assignment. Takes constant time.
        template<int SizeMin>
        void operator=(DynamicList<T, SizeMin>&& lst);

        //- Move assignment. Takes constant time.
        void operator=(SortableList<T>&& lst);


    // Istream operator

        //- Read List from Istream, discarding contents of existing List
        friend Istream& operator>> <T>
        (
            Istream& is,
            List<T>& L
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "List.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
