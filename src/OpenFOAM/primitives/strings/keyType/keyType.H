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
    Foam::keyType

Description
    A class for handling keywords in dictionaries.

    A keyType is the keyword of a dictionary.
    It differs from word in that it also accepts patterns (regular expressions).
    It is very similar to wordRe, but doesn't store the regular expression
    separately.

SourceFiles
    keyType.C

\*---------------------------------------------------------------------------*/

#ifndef keyType_H
#define keyType_H

#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class keyType;
Istream& operator>>(Istream& is, keyType& kw);
Ostream& operator<<(Ostream& os, const keyType& kw);


/*---------------------------------------------------------------------------*\
                           Class keyType Declaration
\*---------------------------------------------------------------------------*/

class keyType
:
    public word
{
    // Private data

        //- Is the keyType a pattern (regular expression)
        bool isPattern_;

    // Private Member Functions

        //- No assignment where we cannot determine string/word type
        void operator=(const std::string&) = delete;

public:

    // Static data members

        //- An empty keyType
        static const keyType null;


    // Constructors

        //- Construct null
        inline keyType();

        //- Copy construct, retaining type (literal or regex)
        inline keyType(const keyType& s);

        //- Copy construct from word. Not treated as a regular expression
        inline keyType(const word& s);

        //- Copy construct from string. Treat as regular expression.
        inline keyType(const string& s);

        //- Construct as copy of character array.
        //  Not treated as a regular expression
        inline keyType(const char* s);

        //- Copy construct from std::string with specified treatment
        inline keyType(const std::string& s, const bool isPattern);

        //- Move construct, retaining type (literal or regex)
        inline keyType(keyType&& s);

        //- Move construct from word. Not treated as a regular expression
        inline keyType(word&& s);

        //- Move construct from string. Treat as regular expression.
        inline keyType(string&& s);

        //- Move construct from std::string with specified treatment
        inline keyType(std::string&& s, const bool isPattern);

        //- Construct from Istream
        //  Treat as regular expression if surrounded by quotation marks.
        keyType(Istream& is);


    // Member functions

        //- Is this character valid for a keyType?
        //  This is largely identical with what word accepts, but also
        //  permit brace-brackets, which are valid for some regexs.
        inline static bool valid(char c);

        //- Treat as a pattern rather than a literal string?
        inline bool isPattern() const;

        //- Swap contents
        inline void swap(keyType& s);

        //- Smart match as regular expression or as a string.
        //  Optionally force a literal match only
        bool match(const std::string& text, bool literal = false) const;


    // Member operators

        //- Perform smart match on text
        inline bool operator()(const std::string& text) const;


        //- Copy assignment, retaining type (literal or regex)
        inline void operator=(const keyType& s);

        //- Assign as word, not treated as a regular expression.
        inline void operator=(const word& s);

        //- Assign from Foam::string as regular expression
        inline void operator=(const string& s);

        //- Assign as word, not treated as a regular expression.
        inline void operator=(const char* s);

        //- Move assignment, retaining type (literal or regex)
        inline void operator=(keyType&& s);


    // IOstream operators

        friend Istream& operator>>(Istream& is, keyType& kw);
        friend Ostream& operator<<(Ostream& os, const keyType& kw);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "keyTypeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
