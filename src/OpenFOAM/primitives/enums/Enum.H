/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2019 OpenCFD Ltd.
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
    Enum.C
    EnumI.H

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

/*---------------------------------------------------------------------------*\
                            Class Enum Declaration
\*---------------------------------------------------------------------------*/

template<class EnumType>
class Enum
{
    // Private Member Data

        //- The names for the enum
        List<word> keys_;

        //- The values for the enum
        List<int> vals_;


public:

    //- The type of enumeration represented by the Enum
    typedef EnumType value_type;

    // Normally only enum, but be generous and allow integrals as well
    static_assert
    (
        std::is_enum<EnumType>::value || std::is_integral<EnumType>::value,
        "Enum must be enum or an integral type"
    );


    // Constructors

        //- Construct from a values/names list.
        //  Duplicate values are permitted (eg, for aliases).
        //  Duplicate names are permitted, but won't make much sense.
        explicit Enum
        (
            std::initializer_list<std::pair<EnumType, const char*>> list
        );


    // Member Functions

    // Access

        //- The number of name/value pairs for the enumeration.
        inline label size() const;

        //- The list of enum names, in construction order. Same as toc()
        inline const List<word>& names() const;

        //- The list of enum names, in construction order. Same as names()
        inline const List<word>& toc() const;

        //- The sorted list of enum names.
        List<word> sortedToc() const;

        //- The list of enum values, in construction order.
        inline const List<int>& values() const;


    // Query

        //- Find the index of the given name.
        //  \return position in list or -1 if not found.
        inline label find(const word& enumName) const;

        //- Find the first index of given enumeration.
        //  \return position in list or -1 if not found.
        inline label find(const EnumType e) const;

        //- Test if there is an enumeration corresponding to the given name.
        inline bool found(const word& enumName) const;

        //- Test if there is a name corresponding to the given enumeration.
        inline bool found(const EnumType e) const;

        //- The enumeration corresponding to the given name.
        //  FatalError if not found.
        EnumType get(const word& enumName) const;

        //- The enumeration corresponding to the given name.
        //  \return The enumeration or default if not found.
        EnumType get(const word& enumName, const EnumType deflt) const;

        //- The name corresponding to the given enumeration.
        //  Return an empty word if not found.
        inline const word& get(const EnumType e) const;


    // Read

        //- Get the key in the dictionary and return the corresponding
        //- enumeration element based on its name.
        //  FatalError if anything is incorrect.
        EnumType get
        (
            const word& key,        //!< Lookup key. Uses LITERAL (not REGEX)
            const dictionary& dict  //!< dictionary
        ) const;

        //- Find the key in the dictionary and return the corresponding
        //- enumeration element based on its name.
        //
        //  \return The value found or default if not found in dictionary.
        //  FatalError (or Warning) if the enumerated name was incorrect.
        EnumType getOrDefault
        (
            const word& key,        //!< Lookup key. Uses LITERAL (not REGEX)
            const dictionary& dict, //!< dictionary
            const EnumType deflt,   //!< fallback if not found
            const bool failsafe = false  //!< Warn only on bad enumeration
        ) const;

        //- Find entry and assign to T val.
        //- FatalIOError if it is found and the number of tokens is incorrect,
        //- or it is mandatory and not found.
        //
        //  \return true if the entry was found.
        bool readEntry
        (
            const word& key,        //!< Lookup key. Uses LITERAL (not REGEX)
            const dictionary& dict, //!< dictionary
            EnumType& val,          //!< the value to read into
            bool mandatory = true   //!< the keyword is mandatory
        ) const;

        //- Find an entry if present, and assign to T val.
        //- FatalIOError if it is found and the number of tokens is incorrect.
        //  Default search: non-recursive with patterns.
        //
        //  \return true if the entry was found.
        bool readIfPresent
        (
            const word& key,        //!< Lookup key. Uses LITERAL (not REGEX)
            const dictionary& dict, //!< dictionary
            EnumType& val           //!< the value to read into
        ) const;


    // IO

        //- Read a word from Istream and return the corresponding enumeration
        EnumType read(Istream& is) const;

        //- Write the name representation of the enumeration to an Ostream
        //  A noop if the enumeration wasn't found.
        inline void write(const EnumType e, Ostream& os) const;

        //- Write the names as a list to an Ostream.
        //  Default is without line-breaks.
        inline Ostream& writeList(Ostream& os, const label shortLen=0) const;


    // Member Operators

        //- Return the enumeration corresponding to the given name
        //  FatalError if the name is not found.
        //  Identical to single parameter get()
        inline EnumType operator[](const word& enumName) const;

        //- Return the first name corresponding to the given enumeration.
        //  Returns an empty word on failure.
        //  Identical to single parameter get()
        inline const word& operator[](const EnumType e) const;

        //- Return the enumeration corresponding to the given name or
        //- deflt if the name is not found.
        inline EnumType operator()
        (
            const word& enumName,
            const EnumType deflt
        ) const;


    // Housekeeping

        //- Find the key in the dictionary and return the corresponding
        //- enumeration element based on its name.
        //
        //  \return The value found or default if not found in dictionary.
        //  FatalError (or Warning) if the enumerated name was incorrect.
        EnumType lookupOrDefault
        (
            const word& key,        //!< Lookup key. Uses LITERAL (not REGEX)
            const dictionary& dict, //!< dictionary
            const EnumType deflt,   //!< fallback if not found
            const bool failsafe = false //!< Warn only on bad enumeration
        ) const
        {
            return getOrDefault(key, dict, deflt, failsafe);
        }

        //- Deprecated(2018-10) same as two-parameter get()
        //  \deprecated(2018-10) - use two-parameter get() method
        EnumType FOAM_DEPRECATED_FOR(2018-10, "get() method")
        lookup(const word& key, const dictionary& dict) const
        {
            return get(key, dict);
        }

        //- Deprecated(2018-10) lookupOrDefault with warnings instead of error.
        //  \deprecated(2018-10) - use getOrDefault() with failsafe option
        EnumType FOAM_DEPRECATED_FOR(2018-10, "getOrDefault() method")
        lookupOrFailsafe
        (
            const word& key,
            const dictionary& dict,
            const EnumType deflt
        ) const
        {
            return getOrDefault(key, dict, deflt, true);
        }
};


// Ostream Operator

//- Write enumeration names, without line-breaks (ie, FlatOutput)
template<class EnumType>
inline Ostream& operator<<(Ostream& os, const Enum<EnumType>& list);


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
