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

Class
    Foam::Detail::HashTablePair

Description
    Internal storage type for HashTable

SourceFiles
    HashTableDetail.H

\*---------------------------------------------------------------------------*/

#ifndef HashTableDetail_H
#define HashTableDetail_H

#include "zero.H"
#include <utility>
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class Ostream;

namespace Detail
{

/*---------------------------------------------------------------------------*\
                    Class Detail::HashTablePair Declaration
\*---------------------------------------------------------------------------*/

//- Internal storage type for HashTable entries
//  Structure with a (K,V) tuple and a linked-list for collisions
//  Could store key/val as std::pair, but no particular advantage
//  unless the iterator dereference type changes.
template<class T, class Key>
struct HashTablePair
{
    // Types

        //- Type of key
        typedef Key key_type;

        //- Type of content
        typedef T mapped_type;

        //- This struct stores a value
        static constexpr bool stores_value()
        {
            return true;
        }


    // Member Data

        //- The lookup key
        key_type key_;

        //- The mapped data
        mapped_type val_;

        //- Addressing (next in collision list)
        HashTablePair* next_;


    // Member Functions

        //- No copy construct
        HashTablePair(const HashTablePair&) = delete;

        //- No copy assignment
        void operator=(const HashTablePair&) = delete;


        //- Construct from key, value, next pointer
        HashTablePair
        (
            const key_type& key,
            const mapped_type& val,
            HashTablePair* next
        )
        :
            key_(key),
            val_(val),
            next_(next)
        {}


        //- The key
        const key_type& key() const
        {
            return key_;
        }

        //- The mapped value
        const mapped_type& mapped() const
        {
            return val_;
        }

        //- The mapped value
        mapped_type& mapped()
        {
            return val_;
        }

        //- Write (key, val) pair - for pointer types
        template<class TypeT = T>
        typename std::enable_if<std::is_pointer<TypeT>::value, void>::type
        print(Ostream& os) const
        {
            os << key_;

            if (val_)
            {
                os << ' ' << *val_;
            }
        }

        //- Write (key, val) pair - for non-pointer types
        template<class TypeT = T>
        typename std::enable_if<!std::is_pointer<TypeT>::value, void>::type
        print(Ostream& os) const
        {
            os << key_ << ' ' << val_;
        }
};


/*---------------------------------------------------------------------------*\
                   Class Detail::HashTableSingle Declaration
\*---------------------------------------------------------------------------*/

//- Internal storage type for HashSet entries
//  Structure with a single (K) value and a linked-list for collisions
template<class T, class Key>
struct HashTableSingle
{
    // Types

        //- Type of key
        typedef Key key_type;

        //- Type of content
        typedef Foam::zero::null mapped_type;

        //- Content storage type to the entry
        typedef key_type value_type;

        //- This struct does not store a value
        static constexpr bool stores_value()
        {
            return false;
        }

    // Member Data

        //- The lookup key == content
        key_type key_;

        //- Addressing (next in collision list)
        HashTableSingle* next_;


    // Member Functions

        //- No copy construct
        HashTableSingle(const HashTableSingle&) = delete;

        //- No copy assignment
        void operator=(const HashTableSingle&) = delete;


        //- Construct from key, (ununsed) value, next pointer
        HashTableSingle
        (
            const key_type& key,
            const mapped_type&,
            HashTableSingle* next
        )
        :
            key_(key),
            next_(next)
        {}


        //- The key
        const key_type& key() const
        {
            return key_;
        }

        //- The mapped value (dummy)
        const mapped_type& mapped() const
        {
            return Foam::zero::null::dummy;
        }

        //- The mapped value (dummy)
        mapped_type& mapped()
        {
            return Foam::zero::null::dummy;
        }

        //- Write the key. There is no value
        void print(Ostream& os) const
        {
            os << key_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Detail
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
