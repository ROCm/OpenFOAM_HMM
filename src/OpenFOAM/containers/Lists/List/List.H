/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
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
#include "one.H"
#include "SLListFwd.H"

#include <initializer_list>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class Istream;
class Ostream;

template<class T> class List;
template<class T, unsigned N> class FixedList;
template<class T, int SizeMin> class DynamicList;

template<class T> class PtrList;
template<class T> class SortableList;
template<class T, class Addr> class IndirectListBase;

template<class T> Istream& operator>>(Istream& is, List<T>& list);

// Common list types
typedef List<bool> boolList;        //!< A List of bools
typedef List<char> charList;        //!< A List of chars
typedef List<label> labelList;      //!< A List of labels


/*---------------------------------------------------------------------------*\
                           Class List Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class List
:
    public UList<T>
{
    // Private Member Functions

        //- Allocate list storage
        inline void doAlloc();

        //- Reallocate list storage to the given size
        inline void reAlloc(const label len);

        //- Copy list of given type.
        template<class List2>
        inline void copyList(const List2& list);

        //- Change allocation size of List. Backend for resize/setSize.
        void doResize(const label newSize);

        //- Construct given begin/end iterators and number of elements
        //  Since the size is provided, the end iterator is actually ignored.
        template<class InputIterator>
        inline List
        (
            InputIterator begIter,
            InputIterator endIter,
            const label len
        );


public:

    // Related types

        //- Declare type of subList
        typedef SubList<T> subList;


    // Static Member Functions

        //- Return a null List
        inline static const List<T>& null();


    // Constructors

        //- Null constructor
        inline constexpr List() noexcept;

        //- Construct with given size
        explicit List(const label len);

        //- Construct with given size and value for all elements
        List(const label len, const T& val);

        //- Construct with given size initializing all elements to zero
        List(const label len, const zero);

        //- Construct with length=1, copying the value as the only content
        List(const one, const T& val);

        //- Construct with length=1, moving the value as the only content
        List(const one, T&& val);

        //- Construct with length=1, initializing content to zero
        List(const one, const zero);

        //- Copy construct from list
        List(const List<T>& a);

        //- Copy construct contents from list
        explicit List(const UList<T>& a);

        //- Construct as copy or re-use as specified
        List(List<T>& a, bool reuse);

        //- Construct as subset
        List(const UList<T>& list, const labelUList& mapAddressing);

        //- Construct given begin/end iterators.
        //  Uses std::distance for the size.
        template<class InputIterator>
        List(InputIterator begIter, InputIterator endIter);

        //- Construct as copy of FixedList\<T, N\>
        template<unsigned N>
        explicit List(const FixedList<T, N>& list);

        //- Construct as copy of PtrList<T>
        explicit List(const PtrList<T>& list);

        //- Construct as copy of SLList<T>
        explicit List(const SLList<T>& list);

        //- Construct as copy of IndirectList contents
        template<class Addr>
        explicit List(const IndirectListBase<T, Addr>& list);

        //- Construct from an initializer list
        List(std::initializer_list<T> list);

        //- Move construct from List
        List(List<T>&& list);

        //- Move construct from DynamicList
        template<int SizeMin>
        List(DynamicList<T, SizeMin>&& list);

        //- Move construct from SortableList
        List(SortableList<T>&& list);

        //- Move construct from SLList
        List(SLList<T>&& list);

        //- Construct from Istream
        List(Istream& is);

        //- Clone
        inline autoPtr<List<T>> clone() const;


    //- Destructor
    ~List();


    // Member Functions

        // Edit

            //- Adjust allocated size of list.
            //  The boolList version fills new memory with false.
            inline void resize(const label newSize);

            //- Adjust allocated size of list and set val for new elements
            void resize(const label newSize, const T& val);

            //- Alias for resize(const label)
            inline void setSize(const label newSize);

            //- Alias for resize(const label, const T&)
            inline void setSize(const label newSize, const T& val);

            //- Clear the list, i.e. set size to zero
            inline void clear();

            //- Append an element at the end of the list
            inline void append(const T& val);

            //- Move append an element at the end of the list
            inline void append(T&& val);

            //- Append a List to the end of this list
            inline void append(const UList<T>& list);

            //- Append IndirectList contents at the end of this list
            template<class Addr>
            inline void append(const IndirectListBase<T, Addr>& list);

            //- Transfer the contents of the argument List into this list
            //- and annul the argument list
            void transfer(List<T>& list);

            //- Transfer the contents of the argument List into this list
            //- and annul the argument list
            template<int SizeMin>
            void transfer(DynamicList<T, SizeMin>& list);

            //- Transfer the contents of the argument List into this list
            //- and annul the argument list
            void transfer(SortableList<T>& list);

            //- Return subscript-checked element of UList and resizing the list
            //- if required.
            inline T& newElmt(const label i);


    // Member Operators

        //- Assignment to UList operator. Takes linear time
        void operator=(const UList<T>& a);

        //- Assignment operator. Takes linear time
        void operator=(const List<T>& list);

        //- Assignment to SLList operator. Takes linear time
        void operator=(const SLList<T>& list);

        //- Assignment from IndirectList. Takes linear time
        template<class Addr>
        void operator=(const IndirectListBase<T, Addr>& list);

        //- Assignment to an initializer list
        void operator=(std::initializer_list<T> list);

        //- Assignment of all entries to the given value
        inline void operator=(const T& val);

        //- Assignment of all entries to zero
        inline void operator=(const zero);

        //- Move assignment. Takes constant time
        void operator=(List<T>&& list);

        //- Move assignment. Takes constant time.
        template<int SizeMin>
        void operator=(DynamicList<T, SizeMin>&& list);

        //- Move assignment. Takes constant time.
        void operator=(SortableList<T>&& list);

        //- Move assignment. Takes constant time
        void operator=(SLList<T>&& list);


    // Istream Operator

        //- Read List from Istream, discarding contents of existing List
        friend Istream& operator>> <T>
        (
            Istream& is,
            List<T>& list
        );


    // Housekeeping

        //- No shallowCopy permitted
        void shallowCopy(const UList<T>&) = delete;


    // Special Methods

        //- A bitSet::set() method for a list of bool
        //  Increases size when setting an out-of-bounds value.
        //
        //  \return True if value changed.
        template<class TypeT = T>
        typename std::enable_if<std::is_same<bool, TypeT>::value, bool>::type
        inline set(const label i, bool val = true)
        {
            if (i < 0)
            {
                return false;               // Out-of-bounds: ignore
            }
            else if (i >= this->size())
            {
                if (!val)                   // Unset out-of-bounds: ignore
                {
                    return false;
                }
                this->resize(i+1, false);   // Adjust size for assign, fill 0
            }

            (*this)[i] = val;
            return true;
        }
};


// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

//- Create identity map of the given length with (map[i] == i)
//  Optionally with an alternative start index, so that (map[i] == i+start)
labelList identity(const label len, label start=0);


//- Hashing for List data, which uses Hasher for contiguous data and
//- element-wise incrementally hashing otherwise.
template<class T>
struct Hash<List<T>>
{
    inline unsigned operator()(const UList<T>& obj, unsigned seed=0) const
    {
        if (is_contiguous<T>::value)
        {
            return Hasher(obj.cdata(), obj.size()*sizeof(T), seed);
        }

        for (const T& val : obj)
        {
            seed = Hash<T>()(val, seed);
        }
        return seed;
    }
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
