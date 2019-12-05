/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2018 OpenCFD Ltd.
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
    Foam::HashPtrTable

Description
    A HashTable of pointers to objects of type \<T\>.

SourceFiles
    HashPtrTable.C
    HashPtrTableIO.C

\*---------------------------------------------------------------------------*/

#ifndef HashPtrTable_H
#define HashPtrTable_H

#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations

class Istream;
class Ostream;

template<class T> class autoPtr;
template<class T, class Key, class Hash> class HashPtrTable;

template<class T, class Key, class Hash>
Istream& operator>>(Istream& is, HashPtrTable<T, Key, Hash>& tbl);


/*---------------------------------------------------------------------------*\
                        Class HashPtrTable Declaration
\*---------------------------------------------------------------------------*/

template<class T, class Key=word, class Hash=string::hash>
class HashPtrTable
:
    public HashTable<T*, Key, Hash>
{
    // Private Member Functions

        //- Read from Istream using Istream constructor class
        template<class INew>
        void readIstream(Istream& is, const INew& inew);

        //- Read from dictionary using Istream constructor class
        template<class INew>
        void read(const dictionary& dict, const INew& inew);


public:

    //- The template instance used for this table
    typedef HashPtrTable<T, Key, Hash> this_type;

    //- The template instance used for the parent HashTable
    typedef HashTable<T*, Key, Hash> parent_type;

    using iterator = typename parent_type::iterator;
    using const_iterator = typename parent_type::const_iterator;


    // Constructors

        //- Construct null with default table capacity
        inline HashPtrTable();

        //- Construct given initial table capacity
        inline explicit HashPtrTable(const label size);

        //- Construct from Istream using given Istream constructor class
        template<class INew>
        HashPtrTable(Istream& is, const INew& inew);

        //- Construct from Istream using default Istream constructor class
        HashPtrTable(Istream& is);

        //- Construct from dictionary with default dictionary constructor class
        explicit HashPtrTable(const dictionary& dict);

        //- Copy construct
        HashPtrTable(const this_type& ht);

        //- Move construct
        HashPtrTable(this_type&& ht);


    //- Destructor
    ~HashPtrTable();


    // Member Functions

    // Edit

        //- Remove entry specified by given iterator.
        //  Includes a safeguard against the end-iterator.
        //
        //  \return entry as an encapsulated pointer.
        autoPtr<T> remove(iterator& iter);

        //- Remove entry specified by given key.
        //
        //  \return entry as an encapsulated pointer.
        autoPtr<T> remove(const Key& key);

        //- Erase entry specified by given iterator and delete the
        //- allocated pointer.
        //  Includes a safeguard against the end-iterator.
        //
        //  \return True if item was removed
        bool erase(iterator& iter);

        //- Erase entry specified by given key and delete the
        //- allocated pointer.
        //
        //  \return True if item was removed
        bool erase(const Key& key);

        //- Clear all entries from table and delete any allocated pointers
        void clear();

        //- Write
        void write(Ostream& os) const;


    // Member Operators

        //- Copy assignment
        void operator=(const this_type& rhs);

        //- Move assignment
        void operator=(this_type&& rhs);


    // IOstream Operators

        friend Istream& operator>> <T, Key, Hash>
        (
            Istream& is,
            HashPtrTable<T, Key, Hash>& tbl
        );


    // Housekeeping

        //- No insert() with raw pointers (potential memory leaks).
        //- Use insert() with autoPtr or set()
        inline bool insert(const Key&, T*) = delete;

        //- Insert a new entry, not overwriting existing entries.
        //
        //  \return True if the entry inserted (not previously in table)
        inline bool insert(const Key& key, autoPtr<T>& aptr);

        //- Insert a new entry, not overwriting existing entries.
        //
        //  \return True if the entry inserted (not previously in table)
        inline bool insert(const Key& key, autoPtr<T>&& aptr);

        //- Assign a new entry, overwriting existing entries.
        inline bool set(const Key& key, T* ptr);

        //- Assign a new entry, overwriting existing entries.
        inline bool set(const Key& key, autoPtr<T>& aptr);

        //- Assign a new entry, overwriting existing entries.
        inline bool set(const Key& key, autoPtr<T>&& aptr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "HashPtrTableI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HashPtrTable.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
