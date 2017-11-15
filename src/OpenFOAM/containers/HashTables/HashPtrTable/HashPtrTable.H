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

class Istream;
class Ostream;

// Forward declaration of friend functions and operators

template<class T, class Key, class Hash> class HashPtrTable;

template<class T, class Key, class Hash>
Istream& operator>>(Istream& is, HashPtrTable<T, Key, Hash>& tbl);

template<class T, class Key, class Hash>
Ostream& operator<<(Ostream& os, const HashPtrTable<T, Key, Hash>& tbl);


/*---------------------------------------------------------------------------*\
                        Class HashPtrTable Declaration
\*---------------------------------------------------------------------------*/

template<class T, class Key=word, class Hash=string::hash>
class HashPtrTable
:
    public HashTable<T*, Key, Hash>
{
    // Private Member Functions

        //- Read from Istream using given Istream constructor class
        template<class INew>
        void read(Istream& is, const INew& inewt);

        //- Read from dictionary using given dictionary constructor class
        template<class INew>
        void read(const dictionary& dict, const INew& inewt);


public:

    //- The template instance used for this table
    typedef HashPtrTable<T, Key, Hash> this_type;

    //- The template instance used for the parent HashTable
    typedef HashTable<T*, Key, Hash> parent_type;

    using iterator = typename parent_type::iterator;
    using const_iterator = typename parent_type::const_iterator;


    // Constructors

        //- Construct null with default table capacity
        HashPtrTable();

        //- Construct given initial table capacity
        explicit HashPtrTable(const label size);

        //- Construct from Istream using given Istream constructor class
        template<class INew>
        HashPtrTable(Istream& is, const INew& inewt);

        //- Construct from Istream using default Istream constructor class
        HashPtrTable(Istream& is);

        //- Construct from dictionary with default dictionary constructor class
        explicit HashPtrTable(const dictionary& dict);

        //- Construct as copy
        HashPtrTable(const this_type& ht);


    //- Destructor
    ~HashPtrTable();


    // Member Functions

      // Edit

        //- Remove and return the pointer specified by given iterator.
        //  Includes a safeguard against the end-iterator.
        T* remove(iterator& iter);

        //- Erase an entry specified by given iterator.
        //  Includes a safeguard against the end-iterator.
        bool erase(iterator& iter);

        //- Erase an entry specified by the given key
        bool erase(const Key& key);

        //- Clear all entries from table and deleting any allocated pointers
        void clear();

        //- Write
        void write(Ostream& os) const;


    // Member Operators

        //- Copy assignment
        void operator=(const this_type& rhs);


    // IOstream Operators

        friend Istream& operator>> <T, Key, Hash>
        (
            Istream& is,
            HashPtrTable<T, Key, Hash>& tbl
        );

        friend Ostream& operator<< <T, Key, Hash>
        (
            Ostream& os,
            const HashPtrTable<T, Key, Hash>& tbl
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "HashPtrTable.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
