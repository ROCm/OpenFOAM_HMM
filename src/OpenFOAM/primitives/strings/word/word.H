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
    A class for handling words, derived from Foam::string.

    A word is a string of characters without whitespace, quotes, slashes,
    semicolons or brace brackets. Words are delimited by whitespace.

SourceFiles
    wordI.H
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
        //  Trips an abort on invalid characters for debug 2 or greater
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

        //- Copy construct
        inline word(const word& w);

        //- Construct as copy of character array
        inline word(const char* s, const bool doStripInvalid=true);

        //- Construct as copy with a maximum number of characters
        inline word
        (
            const char* s,
            const size_type len,
            const bool doStripInvalid
        );

        //- Construct as copy of Foam::string
        inline word(const string& s, const bool doStripInvalid=true);

        //- Construct as copy of std::string
        inline word(const std::string& s, const bool doStripInvalid=true);

        //- Move construct
        inline word(word&& w);

        //- Move construct from Foam::string
        inline word(string&& s, const bool doStripInvalid=true);

        //- Move construct from std::string
        inline word(std::string&& s, const bool doStripInvalid=true);

        //- Construct from Istream
        word(Istream& is);


    // Member functions

        //- Is this character valid for a word?
        inline static bool valid(char c);

        //- Construct validated word (no invalid characters).
        //  Optionally prefix any leading digit with '_' to have words
        //  that work nicely as dictionary keywords.
        static word validate(const std::string& s, const bool prefix=false);


      // File-like functions

        //- Return word without extension (part before last .)
        word lessExt() const;

        //- Return file name extension (part after last .)
        word ext() const;

        //- Append a '.' and the ending, and return the object.
        //  The '.' and ending will not be added when the ending is empty,
        //  or when the file name is empty or ended with a '/'.
        word& ext(const word& ending);

        //- Return true if it has an extension or simply ends with a '.'
        inline bool hasExt() const;

        //- Return true if the extension is the same as the given ending.
        bool hasExt(const word& ending) const;

        //- Return true if the extension matches the given ending.
        bool hasExt(const wordRe& ending) const;

        //- Remove extension, returning true if string changed.
        inline bool removeExt();


    // Member operators

      // Assignment

        //- Copy assignment, no character validation required
        inline void operator=(const word& w);

        //- Copy assignment from Foam::string, stripping invalid characters
        inline void operator=(const string& s);

        //- Copy assignment from std::string, stripping invalid characters
        inline void operator=(const std::string& s);

        //- Copy, stripping invalid characters
        inline void operator=(const char* s);

        //- Move assignment
        inline void operator=(word&& w);

        //- Move assignment from Foam::string, stripping invalid characters
        inline void operator=(string&& s);

        //- Move assignment from std::string, stripping invalid characters
        inline void operator=(std::string&& s);


    // IOstream operators

        friend Istream& operator>>(Istream& is, word& w);
        friend Ostream& operator<<(Ostream& os, const word& w);
};


// Global Operators

//- Join words as camelCase, capitalizing the first letter of b.
//  No effect if either argument is empty.
word operator&(const word& a, const word& b);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wordI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
