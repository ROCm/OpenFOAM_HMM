/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::Switch

Description
    A simple wrapper around bool so that it can be read as a word:
    on/off, true/false, yes/no or y/n.

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

Istream& operator>>(Istream&, Switch&);
Ostream& operator<<(Ostream&, const Switch&);

class dictionary;

/*---------------------------------------------------------------------------*\
                           Class Switch Declaration
\*---------------------------------------------------------------------------*/

class Switch
{
    // Private data

        bool bool_;

        word word_;

public:
    // Private member functions

        //- Return a bool representation of a word
        static bool asBool(const word&);

        //- Return a word representation of a bool value
        static word asWord(const bool);

    // Constructors

        //- Construct from bool
        Switch(const bool value)
        :
            bool_(value),
            word_(asWord(value))
        {}

        //- Construct from word
        Switch(const word& value)
        :
            bool_(asBool(value)),
            word_(value)
        {}

        //- Construct from Istream
        Switch(Istream& is);

        //- Construct from dictionary, supplying default value so that if the
        //  value is not found, it is added into the dictionary.
        static Switch lookupOrAddToDict
        (
            const word&,
            dictionary&,
            const Switch& defaultValue = false
        );


    // Member Operators

        Switch& operator=(const Switch& rhs)
        {
            bool_ = rhs.bool_;
            word_ = rhs.word_;

            return *this;
        }

        operator bool() const
        {
            return bool_;
        }

        operator const word&() const
        {
            return word_;
        }


    // Member fuctions

        //- Update the value of the Switch if it is found in the dictionary
        bool readIfPresent(const word&, const dictionary&);


    // IOstream Operators

        friend Istream& operator>>(Istream&, Switch&);
        friend Ostream& operator<<(Ostream&, const Switch&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
