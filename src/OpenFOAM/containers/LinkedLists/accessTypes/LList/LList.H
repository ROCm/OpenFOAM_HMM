/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
    Foam::LList

Description
    Template class for non-intrusive linked lists.

SourceFiles
    LList.C
    LListIO.C

\*---------------------------------------------------------------------------*/

#ifndef LList_H
#define LList_H

#include "label.H"
#include <initializer_list>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

class Istream;
class Ostream;

template<class LListBase, class T> class LList;

template<class LListBase, class T>
Istream& operator>>
(
    Istream& is,
    LList<LListBase, T>& lst
);

template<class LListBase, class T>
Ostream& operator<<
(
    Ostream& os,
    const LList<LListBase, T>&  lst
);


/*---------------------------------------------------------------------------*\
                            Class LList Declaration
\*---------------------------------------------------------------------------*/

template<class LListBase, class T>
class LList
:
    public LListBase
{
public:

    // STL type definitions

        //- Type of values stored.
        typedef T value_type;

        //- Pointer for value_type
        typedef T* pointer;

        //- Const pointer for value_type
        typedef const T* const_pointer;

        //- Reference for value_type
        typedef T& reference;

        //- Const reference for value_type
        typedef const T& const_reference;

        //- The type that can represent the container size
        typedef label size_type;

        //- The difference between iterator objects
        typedef label difference_type;


    // Forward declaration of STL iterators

        class iterator;
        class const_iterator;

        using base_iterator = typename LListBase::iterator;
        using const_base_iterator = typename LListBase::const_iterator;


    //- The storage of T with linked nodes
    struct link
    :
        public LListBase::link
    {
        //- Stored object
        T obj_;

        //- Copy construct from given object
        link(const T& obj)
        :
            obj_(obj)
        {}

        //- Move construct from given object
        link(T&& obj)
        :
            obj_(std::move(obj))
        {}


        //- Dereference LListBase::link to obtain address of stored object
        static constexpr T* ptr(typename LListBase::link* node)
        {
            return &(static_cast<link*>(node)->obj_);
        }

        //- Dereference LListBase::link to obtain address of stored object
        static constexpr const T* ptr(const typename LListBase::link* node)
        {
            return &(static_cast<const link*>(node)->obj_);
        }

        //- Dereference LListBase::link to obtain the stored object
        static constexpr T& ref(typename LListBase::link* node)
        {
            return static_cast<link*>(node)->obj_;
        }

        //- Dereference LListBase::link to obtain the stored object
        static constexpr const T& ref(const typename LListBase::link* node)
        {
            return static_cast<const link*>(node)->obj_;
        }
    };


    // Constructors

        //- Default construct
        LList() = default;

        //- Construct and copy insert the initial T item
        explicit LList(const T& item)
        {
            this->insert(item);
        }

        //- Construct and move insert the initial T item
        explicit LList(T&& item)
        {
            this->insert(std::move(item));
        }

        //- Construct from Istream
        explicit LList(Istream& is);

        //- Copy construct
        LList(const LList<LListBase, T>& lst);

        //- Move construct
        LList(LList<LListBase, T>&& lst);

        //- Copy construct from an initializer list
        LList(std::initializer_list<T> lst);


    //- Destructor
    ~LList();


    // Member Functions

        //- The first entry in the list
        reference first()
        {
            return link::ref(LListBase::first());
        }

        //- The first entry in the list (const access)
        const_reference first() const
        {
            return link::ref(LListBase::first());
        }

        //- The last entry in the list
        reference last()
        {
            return link::ref(LListBase::last());
        }

        //- The last entry in the list (const access)
        const_reference last() const
        {
            return link::ref(LListBase::last());
        }


        //- Add copy at head of list
        void insert(const T& item)
        {
            LListBase::insert(new link(item));
        }

        //- Move construct at head of list
        void insert(T&& item)
        {
            LListBase::insert(new link(std::move(item)));
        }


        //- Add copy at tail of list
        void append(const T& item)
        {
            LListBase::append(new link(item));
        }

        //- Move construct at tail of list
        void append(T&& item)
        {
            LListBase::append(new link(std::move(item)));
        }


        //- Remove and return head
        T removeHead()
        {
            auto p = LListBase::removeHead();
            T obj(std::move(link::ref(p)));
            delete p;
            return obj;
        }

        //- Remove and return element
        T remove(link* item)
        {
            auto p = LListBase::remove(item);
            T obj(std::move(link::ref(p)));
            delete p;
            return obj;
        }

        //- Remove and return element specified by iterator
        T remove(iterator& iter)
        {
            auto p = LListBase::remove(iter);
            T obj(std::move(link::ref(p)));
            delete p;
            return obj;
        }


        //- Delete contents of list
        void clear();

        //- Transfer the contents of the argument into this List
        //  and annul the argument list.
        void transfer(LList<LListBase, T>& lst);


    // Member Operators

        //- Copy assignment
        void operator=(const LList<LListBase, T>& lst);

        //- Move assignment
        void operator=(LList<LListBase, T>&& lst);

        //- Copy assignment from an initializer list
        void operator=(std::initializer_list<T> lst);


    // IOstream Operators

        //- Read list from Istream
        Istream& readList(Istream& is);

        //- Write LList with line-breaks when length exceeds shortLen.
        //  Using '0' suppresses line-breaks entirely.
        Ostream& writeList(Ostream& os, const label shortLen=0) const;

        //- Read list from Istream
        friend Istream& operator>> <LListBase, T>
        (
            Istream&,
            LList<LListBase, T>& lst
        );

        //- Write LList to Ostream with line breaks,
        //- as per writeList with shortLen=-1
        friend Ostream& operator<< <LListBase, T>
        (
            Ostream& os,
            const LList<LListBase, T>& lst
        );


    // STL iterator

        //- An STL-conforming iterator
        class iterator
        :
            public base_iterator
        {
        public:

            //- Construct from base iterator
            iterator(base_iterator iter)
            :
                base_iterator(iter)
            {}

            reference operator*() const
            {
                return link::ref(this->get_node());
            }

            pointer operator->() const
            {
                return link::ptr(this->get_node());
            }

            reference operator()() const
            {
                return operator*();
            }

            iterator& operator++()
            {
                this->next();
                return *this;
            }

            iterator& operator--()
            {
                this->prev();   // May not be implemented
                return *this;
            }
        };


    // STL const_iterator

        //- An STL-conforming const_iterator
        class const_iterator
        :
            public const_base_iterator
        {
        public:

            //- Construct from base iterator
            const_iterator(const_base_iterator iter)
            :
                const_base_iterator(iter)
            {}

            //- Construct from base iterator
            const_iterator(base_iterator iter)
            :
                const_base_iterator(iter)
            {}

            const_reference operator*() const
            {
                return link::ref(this->get_node());
            }

            const_pointer operator->() const
            {
                return link::ptr(this->get_node());
            }

            const_reference operator()() const
            {
                return operator*();
            }

            const_iterator& operator++()
            {
                this->next();
                return *this;
            }

            const_iterator& operator--()
            {
                this->prev();   // May not be implemented
                return *this;
            }
        };


    // STL reverse_iterator

        //- A reverse_iterator, for LListBase classes that support
        //- reverse iteration
        class reverse_iterator
        :
            public base_iterator
        {
        public:

            //- Construct from base iterator
            reverse_iterator(base_iterator iter)
            :
                base_iterator(iter)
            {}

            reference operator*() const
            {
                return link::ref(this->get_node());
            }

            pointer operator->() const
            {
                return link::ptr(this->get_node());
            }

            reverse_iterator& operator++()
            {
                this->prev();   // Only if base iterator is bidirectional
                return *this;
            }

            reverse_iterator& operator--()
            {
                this->next();
                return *this;
            }
        };


    // STL const_reverse_iterator

        //- A const_reverse_iterator, for LListBase classes that support
        //- reverse iteration
        class const_reverse_iterator
        :
            public const_base_iterator
        {
        public:

            //- Construct from base iterator
            const_reverse_iterator(const_base_iterator iter)
            :
                const_base_iterator(iter)
            {}

            const_reference operator*() const
            {
                return link::ref(this->get_node());
            }

            const_pointer operator->() const
            {
                return link::ptr(this->get_node());
            }

            const_reverse_iterator& operator++()
            {
                this->prev();   // Only if base iterator is bidirectional
                return *this;
            }

            const_reverse_iterator& operator--()
            {
                this->next();
                return *this;
            }
        };


        //- Iterator to first item in list with non-const access
        inline iterator begin()
        {
            return LListBase::template iterator_first<base_iterator>();
        }

        //- Iterator to first item in list with const access
        inline const_iterator cbegin() const
        {
            return LListBase::template iterator_first<const_base_iterator>();
        }

        //- Iterator to last item in list with non-const access
        inline reverse_iterator rbegin()
        {
            return LListBase::template iterator_last<base_iterator>();
        }

        //- Iterator to last item in list with const access
        inline const_reverse_iterator crbegin() const
        {
            return LListBase::template iterator_last<const_base_iterator>();
        }

        //- Iterator to first item in list with const access
        inline const_iterator begin() const
        {
            return LListBase::cbegin();
        }

        //- Iterator to last item in list with const access
        inline const_reverse_iterator rbegin() const
        {
            return crbegin();
        }


        //- End of list for forward iterators
        inline const iterator& end()
        {
            return LListBase::template iterator_end<iterator>();
        }

        //- End of list for forward iterators
        inline const const_iterator& cend() const
        {
            return LListBase::template iterator_end<const_iterator>();
        }

        //- End of list for reverse iterators
        inline const reverse_iterator& rend()
        {
            return LListBase::template iterator_rend<reverse_iterator>();
        }

        //- End of list for reverse iterators
        inline const const_reverse_iterator& crend() const
        {
            return LListBase::template iterator_rend<const_reverse_iterator>();
        }

        //- End of list for forward iterators
        inline const const_iterator& end() const
        {
            return cend();
        }

        //- End of list for reverse iterators
        inline const const_reverse_iterator& rend() const
        {
            return crend();
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LList.C"
    #include "LListIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
