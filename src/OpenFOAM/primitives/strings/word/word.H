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
    Foam::word

Description
    A class for handling words, derived from string.

    A word is a string of characters without whitespace, quotes, slashes,
    semicolons or brace brackets. Words are delimited by whitespace.

SourceFiles
    word.C
    wordIO.C

\*---------------------------------------------------------------------------*/

#ifndef word_H
#define word_H

#include "string.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class word;
inline word operator&(const word& a, const word& b);
Istream& operator>>(Istream& is, word& w);
Ostream& operator<<(Ostream& os, const word& w);


/*---------------------------------------------------------------------------*\
                           Class word Declaration
\*---------------------------------------------------------------------------*/

class word
:
    public string
{
    // Private Member Functions

        //- Strip invalid characters from this word
        inline void stripInvalid();


public:

    // Static data members

        static const char* const typeName;
        static int debug;

        //- An empty word
        static const word null;


    // Constructors

        //- Construct null
        inline word();

        //- Construct as copy
        inline word(const word& w);

        //- Construct as copy of character array
        inline word(const char* s, const bool doStripInvalid=true);

        //- Construct as copy with a maximum number of characters
        inline word
        (
            const char* s,
            const size_type,
            const bool doStripInvalid
        );

        //- Construct as copy of string
        inline word(const string& s, const bool doStripInvalid=true);

        //- Construct as copy of std::string
        inline word(const std::string& s, const bool doStripInvalid=true);

        //- Construct from Istream
        word(Istream& is);


    // Member functions

        //- Is this character valid for a word?
        inline static bool valid(char c);

        //- Construct a validated word, in which all invalid characters have
        //  been stripped out and any leading digit is '_'-prefixed.
        static word validated(const std::string& s);


    // Member operators

        // Assignment

            inline void operator=(const word& w);
            inline void operator=(const string& s);
            inline void operator=(const std::string& s);
            inline void operator=(const char* s);


    // Friend Operators

        //- Join word a and b, capitalising the first letter of b
        //  (so-called camelCase)
        friend word operator&(const word& a, const word& b);


    // IOstream operators

        friend Istream& operator>>(Istream& is, word& w);
        friend Ostream& operator<<(Ostream& os, const word& w);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wordI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
