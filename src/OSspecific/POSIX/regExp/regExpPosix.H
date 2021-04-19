/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
    Foam::regExpPosix

Description
    Wrapper around POSIX extended regular expressions.

    The PCRE '(?i)' extension is provided to compile the regular expression
    as being case-insensitive.

SeeAlso
    The manpage regex(7) for more information about POSIX regular expressions.
    These differ somewhat from \c Perl and \c sed regular expressions.

SeeAlso
    Foam::regExp and Foam::regExpCxx

Warning
    This class should not be used directly.
    Use the Foam::regExp typedef instead.

\deprecated
    This class will be superseded by Foam::regExpCxx as compiler support
    for regular expressions continues to improve.

SourceFiles
    regExpPosixI.H
    regExpPosix.C

\*---------------------------------------------------------------------------*/

#ifndef regExpPosix_H
#define regExpPosix_H

#include "regExpCxx.H"
#include <regex.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class StringType> class SubStrings;

/*---------------------------------------------------------------------------*\
                         Class regExpPosix Declaration
\*---------------------------------------------------------------------------*/

class regExpPosix
{
    // Private Data

        //- Compiled regular expression
        regex_t* preg_;

public:

    // Public Types

        //- Type for matches - similar to std::smatch
        typedef SubStrings<std::string> results_type;


    // Static Member Data

        //- Grammar (unused) - for compatibility with Foam::regExpCxx
        static int grammar;


    // Static Member Functions

        //- Test if character is a regex meta-character
        inline static bool is_meta(const char c) noexcept
        {
            return regExpCxx::is_meta(c);
        }

        //- Test if string contains any (unquoted) meta-characters
        inline static bool is_meta
        (
            const std::string& str,
            const char quote = '\\'
        )
        {
            return regExpCxx::is_meta(str, quote);
        }


    // Public Classes

        //- Functor wrapper for testing meta-characters
        using meta = regExpCxx::meta;


    // Constructors

        //- Default construct
        inline regExpPosix() noexcept;

        //- Copy construct - disallowed
        regExpPosix(const regExpPosix&) = delete;

        //- Move construct
        inline regExpPosix(regExpPosix&& rgx) noexcept;

        //- Construct from character array, optionally ignore case
        inline explicit regExpPosix
        (
            const char* pattern,
            const bool ignoreCase = false
        );

        //- Construct from string, optionally ignore case
        inline explicit regExpPosix
        (
            const std::string& pattern,
            const bool ignoreCase = false
        );


    //- Destructor
    inline ~regExpPosix();


    // Member Functions

    // Access

        //- Return true if a precompiled expression does not exist
        inline bool empty() const noexcept;

        //- Return true if a precompiled expression exists
        inline bool exists() const noexcept;

        //- The number of capture groups for a non-empty expression
        inline unsigned ngroups() const;


    // Editing

        //- Clear expression.
        //  \return True if expression had existed prior to the clear.
        bool clear();

        //- Swap contents
        inline void swap(regExpPosix& rgx);

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
        void operator=(const regExpPosix&) = delete;

        //- Move assignment
        inline void operator=(regExpPosix&& rgx);

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

#include "regExpPosixI.H"

#endif

// ************************************************************************* //
