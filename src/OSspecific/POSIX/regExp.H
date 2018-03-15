/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017-2018 OpenCFD Ltd.
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
    Foam::regExp

Description
    Wrapper around POSIX extended regular expressions.

    The PCRE '(?i)' extension is provided to compile the regular expression
    as being case-insensitive.

See also
    The manpage regex(7) for more information about POSIX regular expressions.
    These differ somewhat from \c Perl and \c sed regular expressions.

SourceFiles
    regExpI.H
    regExp.C

\*---------------------------------------------------------------------------*/

#ifndef regExp_H
#define regExp_H

#include <regex.h>
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
template<class String> class SubStrings;


/*---------------------------------------------------------------------------*\
                           Class regExp Declaration
\*---------------------------------------------------------------------------*/

class regExp
{
    // Private data

        //- Precompiled regular expression
        regex_t* preg_;

public:

    //- Type for matches
    typedef SubStrings<std::string> results_type;


    // Static Member Functions

        //- Test if character appears to be a regular expression meta-character
        //  \return true if character is one of the following:
        //  - any character: '.' \n
        //  - quantifiers:  '*', '+', '?' \n
        //  - grouping: '(', '|', ')' \n
        //  - range: '[', ']' \n
        //
        // \note The presence of '{', '}' regex bounds is not considered
        inline static bool meta(const char c);


    // Constructors

        //- Construct null
        inline regExp();

        //- Copy construct - disallowed
        regExp(const regExp&) = delete;

        //- Move construct
        inline regExp(regExp&& rgx);

        //- Construct from character array
        inline explicit regExp(const char* pattern);

        //- Construct from string
        inline explicit regExp(const std::string& pattern);

        //- Construct from character array, optionally ignore case
        inline regExp(const char* pattern, bool ignoreCase);

        //- Construct from string, optionally ignore case
        inline regExp(const std::string& pattern, bool ignoreCase);


    //- Destructor
    inline ~regExp();


    // Member functions

    // Access

        //- Return true if a precompiled expression does not exist
        inline bool empty() const;

        //- Return true if a precompiled expression exists
        inline bool exists() const;

        //- The number of capture groups for a non-empty expression
        inline unsigned ngroups() const;

    // Editing

        //- Clear expression.
        //  \return True if expression had existed prior to the clear.
        bool clear();

        //- Swap contents
        inline void swap(regExp& rgx);

        //- Compile pattern into a regular expression, optionally ignore case.
        //  \return True if the pattern was compiled
        bool set(const char* pattern, bool ignoreCase=false);

        //- Compile pattern into a regular expression, optionally ignore case.
        //  \return True if the pattern was compiled
        bool set(const std::string& pattern, bool ignoreCase=false);

    // Matching/Searching

        //- Find position within the text.
        //  \return The index where it begins or string::npos if not found
        std::string::size_type find(const std::string& text) const;

        //- True if the regex matches the entire text.
        //  The begin-of-line (^) and end-of-line ($) anchors are implicit
        bool match(const std::string& text) const;

        //- True if the regex matches the text, set the matches.
        //  The first group starts at index 1 (0 is the entire match).
        //  The begin-of-line (^) and end-of-line ($) anchors are implicit
        bool match(const std::string& text, results_type& matches) const;

        //- Return true if the regex was found within the text
        inline bool search(const std::string& text) const;


    // Member Operators

        //- Perform match on text
        inline bool operator()(const std::string& text) const;

        //- Copy assignment - disallowed
        void operator=(const regExp&) = delete;

        //- Move assignment
        inline void operator=(regExp&& rgx);

        //- Assign and compile pattern from a character array.
        //  Matching is case sensitive.
        inline void operator=(const char* pattern);

        //- Assign and compile pattern from string.
        //  Matching is case sensitive.
        inline void operator=(const std::string& pattern);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "regExpI.H"

#endif

// ************************************************************************* //
