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
    Foam::Switch

Description
    A simple wrapper around bool so that it can be read as a word:
    true/false, on/off, yes/no, y/n, t/f, or none.

SourceFiles
    Switch.C

\*---------------------------------------------------------------------------*/

#ifndef Switch_H
#define Switch_H

#include "bool.H"
#include "word.H"

// Avoid any pre-processor conflicts with enum names
#undef FALSE
#undef TRUE
#undef NO
#undef YES
#undef OFF
#undef ON
#undef NONE

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class Switch;
class dictionary;

// IOstream Operators
Istream& operator>>(Istream& is, Switch& sw);
Ostream& operator<<(Ostream& is, const Switch& sw);

/*---------------------------------------------------------------------------*\
                           Class Switch Declaration
\*---------------------------------------------------------------------------*/

class Switch
{
public:

    //- Switch enumerations corresponding to common text representations.
    //  \note the values specified here are critical for its proper behaviour.
    //  The lower bit is tested for the true/false condition.
    //  The values correspond to an index into the predefined output names
    //  for the c_str() method.
    enum switchType : unsigned char
    {
        FALSE   = 0 /*!< "false" */, TRUE = 1 /*!< "true" */,
        NO      = 2 /*!< "no" */,    YES  = 3 /*!< "yes" */,
        OFF     = 4 /*!< "off" */,   ON   = 5 /*!< "on" */,
        NONE    = 6 /*!< "none" */,  ANY  = 7 /*!< "any" */,
        INVALID = 8 /*!< "invalid" */
    };


private:

    // Private Data

        //- The logic and enumerated text representation stored in a byte
        unsigned char switch_;


    // Static Member Functions

        //- Return enum value for input string
        static switchType parse(const std::string& str, bool allowBad);


public:

    // Constructors

        //- Null constructible as false
        constexpr Switch() noexcept
        :
            switch_(switchType::FALSE)
        {}

        //- Construct from enumerated value
        constexpr Switch(const switchType sw) noexcept
        :
            switch_(sw)
        {}

        //- Construct from bool
        constexpr Switch(const bool b) noexcept
        :
            switch_(b ? switchType::TRUE : switchType::FALSE)
        {}

        //- Construct from int (treat integer as bool value)
        constexpr Switch(const int i) noexcept
        :
            switch_(i ? switchType::TRUE : switchType::FALSE)
        {}

        //- Construct from float with rounding to zero given by
        //- the tolerance (default: 0.5)
        explicit Switch(const float val, const float tol=0.5);

        //- Construct from double with rounding to zero given by
        //- the tolerance (default: 0.5)
        explicit Switch(const double val, const double tol=0.5);

        //- Construct from string - catches bad input.
        explicit Switch(const std::string& str)
        :
            switch_(parse(str, false))
        {}

        //- Construct from character array - catches bad input.
        explicit Switch(const char* str)
        :
            switch_(parse(str, false))
        {}

        //- Construct from string.
        //  Optionally allow bad words, and catch the error elsewhere
        Switch(const std::string& str, bool allowBad)
        :
            switch_(parse(str, allowBad))
        {}

        //- Construct from character array.
        //  Optionally allow bad words, and catch the error elsewhere
        Switch(const char* str, bool allowBad)
        :
            switch_(parse(str, allowBad))
        {}

        //- Construct from dictionary lookup.
        //  FatalError if anything is incorrect.
        Switch
        (
            const word& key,        //!< Lookup key. Uses LITERAL (not REGEX)
            const dictionary& dict  //!< dictionary
        );

        //- Find the key in the dictionary and return the corresponding
        //- switch value, or the default value.
        //  FatalError if anything is incorrect.
        Switch
        (
            const word& key,        //!< Lookup key. Uses LITERAL (not REGEX)
            const dictionary& dict, //!< dictionary
            const Switch deflt      //!< fallback if not found
        );

        //- Construct from Istream
        explicit Switch(Istream& is);


    // Helpers

        //- Construct from dictionary, supplying default value so that if the
        //- value is not found, it is added into the dictionary.
        static Switch getOrAddToDict
        (
            const word& name,           //!< Lookup key. Uses REGEX!
            dictionary& dict,           //!< dictionary
            const Switch deflt = switchType::FALSE //!< default to add
        );


    // Static Member Functions

        //- A string representation of bool as "false" / "true"
        static const char* name(const bool b) noexcept;


    // Member Functions

        //- True if the Switch represents a valid enumeration
        bool valid() const noexcept;

        //- The underlying enumeration value
        switchType type() const noexcept;

        //- A string representation of the Switch value
        const char* c_str() const noexcept;

        //- A string representation of the Switch value
        std::string str() const;

        //- Update the value of the Switch if it is found in the dictionary
        bool readIfPresent
        (
            const word& name,           //!< Lookup key. Uses REGEX!
            const dictionary& dict      //!< dictionary
        );


    // Member Operators

        //- Conversion to bool
        operator bool() const noexcept
        {
            return (switch_ & 0x1);
        }

        //- Assignment from enumerated value
        Switch& operator=(const switchType sw) noexcept
        {
            switch_ = sw;
            return *this;
        }

        //- Assignment from bool
        Switch& operator=(const bool b) noexcept
        {
            switch_ = (b ? Switch::TRUE : Switch::FALSE);
            return *this;
        }


    // Housekeeping

        //- Construct from dictionary, supplying default value so that if the
        //- value is not found, it is added into the dictionary.
        static Switch lookupOrAddToDict
        (
            const word& name,           //!< Lookup key. Uses REGEX!
            dictionary& dict,           //!< dictionary
            const Switch deflt = switchType::FALSE //!< default to add
        )
        {
            return getOrAddToDict(name, dict, deflt);
        }

        //- Deprecated(2018-03) text representation of the Switch value
        //  \deprecated(2018-03) - use c_str() method
        inline const char* asText() const { return c_str(); };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
