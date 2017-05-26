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

Class
    Foam::NamedEnum

Description
    A NamedEnum is a wrapper around a static list of names that represent
    a particular enumeration. Internally it uses a HashTable for quicker
    lookups.

SourceFiles
    NamedEnum.C

\*---------------------------------------------------------------------------*/

#ifndef NamedEnum_H
#define NamedEnum_H

#include "HashTable.H"
#include "stringList.H"
#include "wordList.H"
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class dictionary;

// Forward declaration
template<class Enum, int> class NamedEnum;

/*---------------------------------------------------------------------------*\
                          Class NamedEnum Declaration
\*---------------------------------------------------------------------------*/

template<class Enum, int nEnum>
class NamedEnum
:
    public HashTable<int>
{
    //- The nEnum must be positive (non-zero)
    static_assert(nEnum > 0, "nEnum must be positive (non-zero)");

    //- The type of HashTable used for the lookup.
    typedef HashTable<int> table_type;


    // Private Member Functions

        //- The names as a list of strings
        template<class StringType>
        static List<StringType> getNamesList();

        //- Disallow default bitwise copy construct
        NamedEnum(const NamedEnum&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const NamedEnum&) = delete;


public:

    //- The type of enumeration wrapped by NamedEnum
    typedef Enum enum_type;


    // Static data members

        //- The set of names corresponding to the enumeration Enum
        static const char* names[nEnum];


    // Constructors

        //- Construct from names
        NamedEnum();


    // Member Functions

        //- Read a word from Istream and return the corresponding
        //  enumeration element
        enum_type read(Istream& is) const;

        //- Write the name representation of the enumeration to an Ostream
        void write(const enum_type e, Ostream& os) const;

        //- Lookup the key in the dictionary and return the corresponding
        // enumeration element based on its name.
        // Fatal if anything is incorrect.
        enum_type lookup
        (
            const word& key,
            const dictionary& dict
        ) const;

        //- Find the key in the dictionary and return the corresponding
        //  enumeration element based on its name.
        //  Return the default value if the key was not found in the dictionary.
        //  Fatal if enumerated name was incorrect.
        enum_type lookupOrDefault
        (
            const word& key,
            const dictionary& dict,
            const enum_type deflt
        ) const;

        //- List of enumerations
        static List<enum_type> enums();

        //- The set of names as a list of strings
        static stringList strings();

        //- The set of names as a list of words
        static wordList words();


    // Member Operators

        //- Return the enumeration element corresponding to the given name
        inline const enum_type operator[](const char* name) const
        {
            return enum_type(table_type::operator[](name));
        }

        //- Return the enumeration element corresponding to the given name
        inline const enum_type operator[](const word& name) const
        {
            return enum_type(table_type::operator[](name));
        }

        //- Return the name of the given enumeration element
        inline const char* operator[](const enum_type e) const
        {
            return names[int(e)];
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NamedEnum.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
