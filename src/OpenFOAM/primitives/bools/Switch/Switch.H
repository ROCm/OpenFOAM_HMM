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
    Foam::Switch

Description
    A simple wrapper around bool so that it can be read as a word:
    true/false, on/off, yes/no, y/n, t/f, or none.

SourceFiles
    Switch.C
    SwitchIO.C

\*---------------------------------------------------------------------------*/

#ifndef Switch_H
#define Switch_H

#include "bool.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class Switch;
class dictionary;

Istream& operator>>(Istream& is, Switch& s);
Ostream& operator<<(Ostream& is, const Switch& s);


/*---------------------------------------------------------------------------*\
                           Class Switch Declaration
\*---------------------------------------------------------------------------*/

class Switch
{
    // Private data

        //- The logic and enumerated text representation stored as a single byte
        unsigned char switch_;

public:

    // Public data types

        // Avoid issues with possible pre-processor defines
        #undef FALSE
        #undef TRUE
        #undef NO
        #undef YES
        #undef OFF
        #undef ON
        #undef NONE
        #undef PLACEHOLDER
        #undef INVALID

        //- The various text representations for a switch value.
        //  These also correspond to the entries in names.
        enum switchType
        {
            FALSE   = 0,  TRUE   = 1,
            NO      = 2,  YES    = 3,
            OFF     = 4,  ON     = 5,
            NONE    = 6,  PLACEHOLDER = 7,
            INVALID
        };

    // Static data members

        //- The set of names corresponding to the switchType enumeration
        //  Includes an extra entry for "invalid".
        static const char* names[INVALID+1];


private:

    // Static Member Functions

        //- Return a switchType representation of an input string
        static switchType asEnum
        (
            const std::string& str,
            const bool allowInvalid
        );


public:

    // Constructors

        //- Construct null as false
        Switch()
        :
            switch_(switchType::FALSE)
        {}

        //- Construct from enumerated value
        Switch(const switchType sw)
        :
            switch_(sw)
        {}

        //- Construct from bool
        Switch(const bool b)
        :
            switch_(b ? switchType::TRUE : switchType::FALSE)
        {}

        //- Construct from integer values (treats integer as bool value)
        Switch(const int i)
        :
            switch_(i ? switchType::TRUE : switchType::FALSE)
        {}

        //- Construct from string - catches bad input.
        explicit Switch(const std::string& str)
        :
            switch_(asEnum(str, false))
        {}

        //- Construct from character array - catches bad input.
        explicit Switch(const char* str)
        :
            switch_(asEnum(std::string(str), false))
        {}

        //- Construct from string.
        //  Optionally allow bad words, and catch the error elsewhere
        Switch(const std::string& str, const bool allowInvalid)
        :
            switch_(asEnum(str, allowInvalid))
        {}

        //- Construct from character array.
        //  Optionally allow bad words, and catch the error elsewhere
        Switch(const char* str, const bool allowInvalid)
        :
            switch_(asEnum(std::string(str), allowInvalid))
        {}

        //- Construct from Istream
        Switch(Istream& is);

        //- Construct from dictionary, supplying default value so that if the
        //  value is not found, it is added into the dictionary.
        static Switch lookupOrAddToDict
        (
            const word& name,
            dictionary& dict,
            const Switch& defaultValue = switchType::FALSE
        );


    // Member Functions

        //- Return true if the Switch has a valid value
        bool valid() const;

        //- Return a text representation of the Switch
        const char* asText() const;

        //- Update the value of the Switch if it is found in the dictionary
        bool readIfPresent(const word& name, const dictionary& dict);


    // Member Operators

        //- Conversion to bool
        operator bool() const
        {
            return (switch_ & 0x1);
        }

        //- Assignment to enumerated value
        void operator=(const switchType sw)
        {
            switch_ = sw;
        }

        //- Assignment to bool
        void operator=(const bool b)
        {
            switch_ = (b ? Switch::TRUE : Switch::FALSE);
        }


    // IOstream Operators

        friend Istream& operator>>(Istream& is, Switch& s);
        friend Ostream& operator<<(Ostream& os, const Switch& s);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
