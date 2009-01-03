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
    Foam::wordRe

Description
    A wordRe is a word, but can also have a regular expression for matching
    words.

Note
    If the string contents are changed - eg, by the operator+=() or by
    string::replace(), etc - it will be necessary to use compile() or
    recompile() to sychronize the regular expression.

    THIS IS STILL A DRAFT -- NOT YET RELEASED FOR GENERAL USE

SourceFiles
    wordRe.C
    wordReIO.C

\*---------------------------------------------------------------------------*/

#ifndef wordRe_H
#define wordRe_H

#include "word.H"
#include "regExp.H"
#include "InfoProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class wordRe;

class Istream;
class Ostream;

Istream& operator>>(Istream&, wordRe&);
Ostream& operator<<(Ostream&, const wordRe&);


/*---------------------------------------------------------------------------*\
                           Class wordRe Declaration
\*---------------------------------------------------------------------------*/

class wordRe
:
    public word
{
    // Private member data

        //- The regular expression
        mutable regExp re_;

public:

        //- Is this a meta character?
        static inline bool meta(char);

        //- Is this character valid for a wordRe
        inline static bool valid(char);

        //- Test string for regular expression meta characters
        static inline bool isPattern(const string&);

    // Constructors

        //- Construct null
        inline wordRe();

        //- Construct as copy
        inline wordRe(const wordRe&);

        //- Construct as copy of word
        inline wordRe(const word&);

        //- Construct as copy of character array
        //  Treat as regular expression specified explicitly.
        inline wordRe(const char*, const bool asPattern=false);

        //- Construct as copy of string.
        //  Treat as regular expression specified explicitly.
        inline wordRe(const string&, const bool asPattern);

        //- Construct as copy of string.
        //  Auto-test for regular expression
        inline wordRe(const string&);

        //- Construct as copy of std::string
        //  Treat as regular expression specified explicitly.
        inline wordRe(const std::string&, const bool asPattern);

        //- Construct as copy of std::string
        //  Auto-test for regular expression
        inline wordRe(const std::string&);

        //- Construct from Istream
        wordRe(Istream&);

    // Destructor

        ~wordRe();


    // Member functions

        //- Should be treated as a match rather than a literal string?
        inline bool isPattern() const;

        //- Create and compile the regular expression
        //  Optionally detect if it appears to be a regular expression
        inline bool compile(const bool detect=false) const;

        //- Recompile an existing regular expression
        inline bool recompile() const;

        //- Frees precompiled regular expression and makes is a literal string.
        //  Optionally strips invalid word characters
        inline void uncompile(const bool doStripInvalid=false) const;

        //- Clear string and precompiled regular expression
        inline void clear();

        //- Smart match as regular expression or as a string
        //  Optionally specify a literal match only
        inline bool match(const string&, bool literalMatch=false) const;

        //- Return a string with quoted meta-characters
        inline string quotemeta() const;

        //- Return info proxy.
        InfoProxy<wordRe> info() const
        {
            return *this;
        }


    // Member operators

        // Assignment

            //- Assign copy
            inline void operator=(const wordRe&);

            //- Copy word, never a regular expression
            inline void operator=(const word&);

            //- Copy string, auto-test for regular expression
            inline void operator=(const string&);

            //- Copy string, auto-test for regular expression
            inline void operator=(const std::string&);

            //- Copy string, auto-test for regular expression
            inline void operator=(const char*);


    // IOstream operators

        friend Istream& operator>>(Istream&, wordRe&);
        friend Ostream& operator<<(Ostream&, const wordRe&);
};


template<>
Ostream& operator<<(Ostream&, const InfoProxy<wordRe>&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wordReI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
