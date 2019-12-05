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

// Forward Declarations
class Istream;
class Ostream;

/*---------------------------------------------------------------------------*\
                           Class wordRe Declaration
\*---------------------------------------------------------------------------*/

class wordRe
:
    public word
{
    // Private Member Data

        //- The regular expression
        regExp re_;


public:

    // Static Data Members

        //- An empty wordRe
        static const wordRe null;


    // Public Data Types

        //- Enumeration with compile options
        //  Note that 'REGEX' is implicit if 'ICASE' is specified alone.
        enum compOption
        {
            LITERAL = 0, //!< String literal
            REGEX   = 1, //!< Regular expression
            ICASE   = 2, //!< Ignore case in regular expression
            NOCASE  = 2, //!< \deprecated(2018-04) Alias for ICASE
            DETECT  = 4, //!< Detect if the string contains meta-characters
            UNKNOWN = 4, //!< Unknown content.
            REGEX_ICASE  = (REGEX|ICASE),   //!< Combined REGEX and ICASE
            DETECT_ICASE = (DETECT|ICASE),  //!< Combined DETECT and ICASE
        };


    // Constructors

        //- Construct null
        inline wordRe();

        //- Copy construct
        inline wordRe(const wordRe& str);

        //- Move construct
        inline wordRe(wordRe&& str);

        //- Construct from keyType, using its compile information
        inline explicit wordRe(const keyType& str);

        //- Copy from character array, treat as a literal
        inline explicit wordRe(const char* str);

        //- Copy from std::string, treat as a literal
        inline explicit wordRe(const std::string& str);

        //- Copy from string, treat as a literal
        inline explicit wordRe(const string& str);

        //- Copy from word, treat as a literal
        inline explicit wordRe(const word& str);

        //- Copy from keyType, use specified compile option
        inline wordRe(const keyType& str, const compOption opt);

        //- Copy from character array, use specified compile option
        inline wordRe(const char* str, const compOption opt);

        //- Copy from std::string, use specified compile option
        inline wordRe(const std::string& str, const compOption opt);

        //- Copy from string, use specified compile option
        inline wordRe(const string& str, const compOption opt);

        //- Copy from word, use specified compile option
        inline wordRe(const word& str, const compOption opt);

        //- Construct from Istream
        //  Words are treated as literals, strings with an auto-test
        explicit wordRe(Istream& is);


    // Member Functions

        //- Is this a meta character?
        inline static bool meta(char c);

        //- Is this character valid for a wordRe?
        //  This is largely identical with what word accepts, but also
        //  permit brace-brackets, which are valid for some regexs.
        inline static bool valid(char c);

        //- Test string for regular expression meta characters
        inline static bool isPattern(const std::string& str);


    // Access

        //- The wordRe is treated as literal string, not as pattern.
        inline bool isLiteral() const;

        //- The wordRe is treated as a pattern, not as literal string.
        inline bool isPattern() const;


    // Infrastructure

        //- Compile the regular expression
        inline bool compile();

        //- Possibly compile the regular expression, with greater control
        inline bool compile(const compOption opt);

        //- Make wordRe a literal again, instead of a regular expression.
        //  Optionally strip invalid word characters.
        inline void uncompile(bool adjust = false);


    // Editing

        //- Copy string, auto-test for regular expression or other options
        inline void set(const std::string& str, const compOption opt = DETECT);

        //- Copy string, auto-test for regular expression or other options
        inline void set(const char* str, const compOption opt = DETECT);

        //- Clear string and regular expression
        inline void clear();

        //- Swap contents. Self-swapping is a no-op
        inline void swap(wordRe& str);


    // Matching/Searching

        //- Smart match as regular expression or as a string.
        //  Optionally force a literal match only
        inline bool match(const std::string& text, bool literal=false) const;


    // Miscellaneous

        //- Return a string with quoted meta-characters
        inline string quotemeta() const;

        //- Output some basic info
        Ostream& info(Ostream& os) const;


    // Member Operators

        //- Perform smart match on text, as per match()
        //  Allows use as a predicate.
        inline bool operator()(const std::string& text) const;


        //- Copy assignment, retaining type (literal or regex)
        //  Self-assignment is a no-op.
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
        //  Self-assignment is a no-op.
        inline void operator=(wordRe&& str);
};


// IOstream Operators

//- Read operator
Istream& operator>>(Istream& is, wordRe& val);

//- Write operator
Ostream& operator<<(Ostream& os, const wordRe& val);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wordReI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
