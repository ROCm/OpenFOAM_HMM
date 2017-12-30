/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::Enum

Description
    Enum is a wrapper around a list of names/values that represent particular
    enumeration values.

SourceFiles
    EnumI.H
    Enum.C

\*---------------------------------------------------------------------------*/

#ifndef Enum_H
#define Enum_H

#include "wordList.H"
#include <initializer_list>
#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declarations
class dictionary;
template<class EnumType> class Enum;

template<class EnumType>
Ostream& operator<<(Ostream& os, const Enum<EnumType>& wrapped);

/*---------------------------------------------------------------------------*\
                            Class Enum Declaration
\*---------------------------------------------------------------------------*/

template<class EnumType>
class Enum
{
    // Private Member Data

        //- The names for the enum
        List<word> names_;

        //- The values for the enum
        List<int> values_;


    // Private Member Functions

        //- The index for the given name. -1 if not found.
        label getIndex(const word& enumName) const;

        //- The index for the given enumeration. -1 if not found.
        label getIndex(const EnumType e) const;

        //- Lookup enumeration corresponding to the given name.
        //  FatalError if not found.
        EnumType getEnum(const word& enumName) const;

        //- Lookup name corresponding to the given enumeration.
        //  Return empty word if not found.
        const word& getName(const EnumType e) const;


        //- Disallow default bitwise copy construct
        Enum(const Enum&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const Enum&) = delete;


public:

    //- The type of enumeration wrapped by Enum
    typedef EnumType value_type;


    // Constructors

        //- Construct from a values/names list.
        //  Duplicate values are permitted (eg, for aliases).
        //  Duplicate names are permitted, but won't make much sense.
        explicit Enum(std::initializer_list<std::pair<EnumType, word>> lst);

        //- Construct from a list of names with values incremented from the
        //- specified start value.
        Enum(const EnumType start, std::initializer_list<word> lst);


    // Member Functions

    // Access

        //- The number of lookup names for the enumeration
        inline label size() const;

        //- The list of enum names, in construction order
        inline const List<word>& names() const;

        //- The list of enum names, in construction order
        inline const List<word>& toc() const;

        //- The sorted list of enum names
        List<word> sortedToc() const;

        //- The list of enum values, in construction order
        inline const List<int>& values() const;


    // Query

        //- Test if there is an enumeration corresponding to the given name.
        inline bool found(const word& enumName) const;

        //- Test if there is a name corresponding to the given enumeration.
        inline bool found(const EnumType e) const;

        //- Test if there is an enumeration corresponding to the given name.
        //  \deprecated Use found() - for compatibility with NamedEnum
        inline bool hasEnum(const word& enumName) const;

        //- Test if there is a name corresponding to the given enumeration.
        //  \deprecated Use found() - for compatibility with NamedEnum
        inline bool hasName(const EnumType e) const;


    // Lookup

        //- Lookup the key in the dictionary and return the corresponding
        //- enumeration element based on its name.
        //  FatalError if anything is incorrect.
        EnumType lookup
        (
            const word& key,
            const dictionary& dict
        ) const;

        //- Find the key in the dictionary and return the corresponding
        //- enumeration element based on its name.
        //  Return the default value if the key was not found in the dictionary.
        //  FatalError if the enumerated name was incorrect.
        EnumType lookupOrDefault
        (
            const word& key,
            const dictionary& dict,
            const EnumType deflt
        ) const;


        //- Find the key in the dictionary and return the corresponding
        //- enumeration element based on its name.
        //  Return the default value if the key was not found in the dictionary
        //  or if the enumerated name was incorrect (emit warning)
        EnumType lookupOrFailsafe
        (
            const word& key,
            const dictionary& dict,
            const EnumType deflt
        ) const;


    // IO

        //- Read a word from Istream and return the corresponding enumeration
        EnumType read(Istream& is) const;

        //- Write the name representation of the enumeration to an Ostream
        //  A noop if the enumeration wasn't found.
        void write(const EnumType e, Ostream& os) const;

        //- Write the names as a list to an Ostream
        inline Ostream& writeList
        (
            Ostream& os,
            const label shortListLen=0
        ) const;


    // Member Operators

        //- Return the enumeration corresponding to the given name
        //  FatalError if the name is not found.
        //  Identical to getEnum()
        inline EnumType operator[](const word& enumName) const;

        //- Return the first name corresponding to the given enumeration.
        //  Returns an empty word on failure.
        //  Identical to getName()
        inline const word& operator[](const EnumType e) const;

        //- Return the enumeration corresponding to the given name or deflt
        //- if the name is not found.
        inline EnumType operator()
        (
            const word& enumName,
            const EnumType deflt
        ) const;


    // IOstream operators

        //- Write names to Ostream, as per writeList() with shortListLen=10
        friend Ostream& operator<< <EnumType>
        (
            Ostream& os,
            const Enum<EnumType>& wrapped
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "EnumI.H"

#ifdef NoRepository
    #include "Enum.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
