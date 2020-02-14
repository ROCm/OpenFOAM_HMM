/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Foam::DLListBase

Description
    Base for doubly-linked lists.

    The iterators associated with the list only have a core functionality
    for navigation, with additional functionality to be added by inheriting
    classes. The node iterators always have a node-pointer as the
    first member data, which allows reinterpret_cast from anything else with
    a nullptr as its first data member.
    The nullObject is such an item (with a nullptr data member).

SourceFiles
    DLListBaseI.H
    DLListBase.C

\*---------------------------------------------------------------------------*/

#ifndef DLListBase_H
#define DLListBase_H

#include "label.H"
#include "uLabel.H"
#include "stdFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class DLListBase Declaration
\*---------------------------------------------------------------------------*/

class DLListBase
{
public:

    //- The structure for a doubly-linked storage node
    struct link
    {
        //- Pointer to prev entry in list
        link* prev_ = nullptr;

        //- Pointer to next entry in list
        link* next_ = nullptr;

        //- Default construct
        link() = default;

        //- Check if the node is registered with the list
        inline bool registered() const;

        //- Deregister the node after removal
        inline void deregister();
    };


private:

    // Private Data

       //- Pointer to first element
       link *first_ = nullptr;

       //- Pointer to last element
       link *last_ = nullptr;

       //- Number of elements in the list
       label size_ = 0;


protected:

    // Protected Member Functions

        //- Factory method to return an iterator end
        //  Simply reinterprets a NullObject as a DLListBase iterator.
        template<class IteratorType>
        inline static const IteratorType& iterator_end();

        //- Factory method to return an iterator reverse end
        //  Simply reinterprets a NullObject as a DLListBase iterator.
        template<class IteratorType>
        inline static const IteratorType& iterator_rend();

        //- Return iterator to first item or end-iterator if list is empty
        //  Removes constness which the caller promises to manage.
        template<class IteratorType>
        inline IteratorType iterator_first() const;

        //- Return iterator to last item or end-iterator if list is empty
        //  Removes constness which the caller promises to manage.
        template<class IteratorType>
        inline IteratorType iterator_last() const;


public:

    // Forward declaration of iterators

        class iterator;
        friend class iterator;

        class const_iterator;
        friend class const_iterator;


    // Generated Methods

        //- Default construct
        DLListBase() = default;

        //- No copy construct
        DLListBase(const DLListBase&) = delete;

        //- No copy assignment
        void operator=(const DLListBase&) = delete;

        //- Destructor
        ~DLListBase() = default;


    // Member Functions

        //- The number of elements in list
        inline label size() const noexcept;

        //- True if the list is empty
        inline bool empty() const noexcept;

        //- Return first entry
        inline link* first();

        //- Return const access to first entry
        inline const link* first() const;

        //- Return last entry
        inline link* last();

        //- Return const access to last entry
        inline const link* last() const;


        //- Add at head of list
        void insert(link* item);

        //- Add at tail of list
        void append(link* item);

        //- Swap this element with the one above unless it is at the top
        bool swapUp(link* item);

        //- Swap this element with the one below unless it is at the bottom
        bool swapDown(link* item);

        //- Remove and return head
        link* removeHead();

        //- Remove and return element
        link* remove(link* item);

        // Remove and return element specified by iterator
        inline link* remove(iterator& iter);

        //- Replace oldLink with newLink and return element
        link* replace(link* oldLink, link* newLink);

        //- Replace oldIter with newItem and return element
        inline link* replace(iterator& oldIter, link* newitem);

        //- Clear the list
        inline void clear();

        //- Swap the contents of the list
        inline void swap(DLListBase& lst);

        //- Transfer the contents of the argument into this list
        //- and annul the argument list.
        inline void transfer(DLListBase& lst);


    // iterator

        //- A primitive non-const node iterator.
        //  Needs to be extended by inheriting classes.
        class iterator
        {
            friend class DLListBase;
            friend class const_iterator;

            //- The selected node.
            //  MUST be the first member for easy comparison between iterators
            //  and for reinterpret_cast from nullObject
            link* node_;

            //- The list being iterated on
            DLListBase* list_;

            //- Copy of the node prev/next pointers (to use after removal)
            link copy_;

        public:

            //- Copy construct
            iterator(const iterator&) = default;

            //- Construct for a node on a list
            inline iterator(DLListBase* list, link* item);

            //- The storage node
            inline link* get_node() const;

            //- Pointing at a valid storage node
            inline bool good() const;

            //- Deprecated(2019-01) Pointing at a valid storage node
            //  \deprecated(2019-01) - use good() method
            FOAM_DEPRECATED_FOR(2019-01, "good() method")
            bool found() const
            {
                return this->good();
            }

            //- Move backward through list
            inline void prev();

            //- Move forward through list
            inline void next();

            inline void operator=(const iterator& iter);

            inline bool operator==(const iterator&) const;
            inline bool operator!=(const iterator&) const;
        };


    // const_iterator

        //- A primitive const node iterator (bidirectional).
        //  Must normally be extended by inheriting classes.
        //  Since this iterator works bidirectionally, it can be used as the
        //  basis for a derived const_reverse_iterator
        class const_iterator
        {
            //- The selected node.
            //  MUST be the first member for easy comparison between iterators
            //  and for reinterpret_cast from nullObject
            const link* node_;

            //- The list being iterated on (as pointer for bitwise copy)
            const DLListBase* list_;

        public:

            //- Copy construct
            const_iterator(const const_iterator&) = default;

            //- Construct for a node on a list
            inline const_iterator(const DLListBase* list, const link* item);

            //- Copy construct from a non-const iterator
            inline const_iterator(const DLListBase::iterator& iter);

            //- The storage node
            inline const link* get_node() const;

            //- Pointing at a valid storage node
            inline bool good() const;

            //- Deprecated(2019-01) Pointing at a valid storage node
            //  \deprecated(2019-01) - use good() method
            FOAM_DEPRECATED_FOR(2019-01, "good() method")
            bool found() const
            {
                return this->good();
            }

            //- Move backward through list
            inline void prev();

            //- Move forward through list
            inline void next();

            //- Copy assignment
            const_iterator& operator=(const const_iterator&) = default;

            inline bool operator==(const const_iterator&) const;
            inline bool operator!=(const const_iterator&) const;
        };


        //- Iterator to first item in list with non-const access
        inline iterator begin();

        //- Iterator to first item in list with const access
        inline const_iterator cbegin() const;

        //- Iterator to last item in list with const access
        //  Note that this is not a const_reverse_iterator, this is the
        //  responsibility of any derived classes.
        inline const_iterator crbegin() const;

        //- End of list for iterators
        inline const iterator& end();

        //- End of list for iterators
        inline const const_iterator& cend() const;

        //- End of list for reverse iterators
        inline const const_iterator& crend() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DLListBaseI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
