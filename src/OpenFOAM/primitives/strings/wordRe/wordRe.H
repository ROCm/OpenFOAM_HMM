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
    Foam::wordRe

Description
    A wordRe is a Foam::word, but can contain a regular expression for
    matching words or strings.

    By default the constructors will generally preserve the argument as a
    string literal and the assignment operators will use the wordRe::DETECT
    compOption to scan the string for regular expression meta characters
    and/or invalid word characters and react accordingly.

    The exceptions are when constructing/assigning from another
    Foam::wordRe (preserve the same type) or from a Foam::word (always
    literal).

Note
    If the string contents are changed - eg, by the operator+=() or by
    string::replace(), etc - it will be necessary to use compile() to
    synchronize the regular expression.

SourceFiles
    wordRe.C

\*---------------------------------------------------------------------------*/

#ifndef wordRe_H
#define wordRe_H

#include "word.H"
#include "regExp.H"
#include "keyType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class wordRe;
class Istream;
class Ostream;

Istream& operator>>(Istream& is, wordRe& w);
Ostream& operator<<(Ostream& os, const wordRe& w);


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

    // Static data members

        //- An empty wordRe
        static const wordRe null;


    // Public data types

        //- Enumeration with compile options
        //  Note that 'REGEX' is implicit if 'NOCASE' is specified alone.
        enum compOption
        {
            LITERAL = 0, //!< Treat as a string literal
            DETECT  = 1, //!< Detect if the string contains meta-characters
            REGEX   = 2, //!< Treat as regular expression
            NOCASE  = 4, //!< Ignore case in regular expression
            DETECT_NOCASE = DETECT|NOCASE,  //!< Combined DETECT and NOCASE
            REGEX_NOCASE  = REGEX|NOCASE    //!< Combined REGEX and NOCASE
        };


        //- Is this a meta character?
        inline static bool meta(const char c);

        //- Is this character valid for a wordRe?
        //  This is largely identical with what word accepts, but also
        //  permit brace-brackets, which are valid for some regexs.
        inline static bool valid(char c);

        //- Test string for regular expression meta characters
        inline static bool isPattern(const std::string& str);


    // Constructors

        //- Construct null
        inline wordRe();

        //- Construct as copy
        inline wordRe(const wordRe& str);

        //- Construct from keyType, using its compile information
        inline explicit wordRe(const keyType& str);

        //- Construct as copy of character array, treat as a literal
        inline explicit wordRe(const char* str);

        //- Construct as copy of std::string, treat as a literal
        inline explicit wordRe(const std::string& str);

        //- Construct as copy of string, treat as a literal
        inline explicit wordRe(const string& str);

        //- Construct as copy of word, treat as a literal
        inline explicit wordRe(const word& str);

        //- Construct from keyType, use specified compile option
        inline wordRe(const keyType& str, const compOption);

        //- Construct as copy of character array, use specified compile option
        inline wordRe(const char* str, const compOption);

        //- Construct as copy of std::string, use specified compile option
        inline wordRe(const std::string& str, const compOption);

        //- Construct as copy of string, use specified compile option
        inline wordRe(const string& str, const compOption);

        //- Construct as copy of word, use specified compile option
        inline wordRe(const word& str, const compOption);

        //- Move construct
        inline wordRe(wordRe&& str);

        //- Construct from Istream
        //  Words are treated as literals, strings with an auto-test
        wordRe(Istream& is);


    // Member functions

      // Access

        //- Treat as a pattern rather than a literal string?
        inline bool isPattern() const;


      // Infrastructure

        //- Compile the regular expression
        inline bool compile() const;

        //- Possibly compile the regular expression, with greater control
        inline bool compile(const compOption) const;

        //- Make wordRe a literal again, instead of a regular expression.
        //  Optionally strip invalid word characters.
        inline void uncompile(const bool doStripInvalid = false) const;


      // Editing

        //- Copy string, auto-test for regular expression or other options
        inline void set(const std::string& str, const compOption = DETECT);

        //- Copy string, auto-test for regular expression or other options
        inline void set(const char* str, const compOption = DETECT);

        //- Clear string and regular expression
        inline void clear();

        //- Swap contents
        inline void swap(wordRe& str);


     // Matching/Searching

        //- Smart match as regular expression or as a string.
        //  Optionally force a literal match only
        inline bool match(const std::string& text, bool literal = false) const;


      // Miscellaneous

        //- Return a string with quoted meta-characters
        inline string quotemeta() const;

        //- Output some basic info
        Ostream& info(Ostream& os) const;


    // Member operators

        //- Perform smart match on text, as per match()
        inline bool operator()(const std::string& text) const;


        //- Copy assignment, retaining type (literal or regex)
        inline void operator=(const wordRe& str);

        //- Copy word, never a regular expression
        inline void operator=(const word& str);

        //- Copy keyType and its type (literal or regex)
        //  Always case sensitive
        inline void operator=(const keyType& str);

        //- Copy string, auto-test for regular expression
        //  Always case sensitive
        inline void operator=(const string& str);

        //- Copy string, auto-test for regular expression
        //  Always case sensitive
        inline void operator=(const std::string& str);

        //- Copy string, auto-test for regular expression
        //  Always case sensitive
        inline void operator=(const char* str);

        //- Move assignment.
        inline void operator=(wordRe&& str);


    // IOstream operators

        friend Istream& operator>>(Istream& is, wordRe& w);
        friend Ostream& operator<<(Ostream& os, const wordRe& w);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wordReI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
