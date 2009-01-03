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
    Foam::regExp

Description
    Wrapper around POSIX extended regular expressions.

    The beginning-of-line (^) and the end-of-line ($) anchors are implicit
    by default.

SeeAlso
    The manpage regex(7) for more information about POSIX regular expressions.
    These differ somewhat from @c Perl and @c sed regular expressions.

SourceFiles
    regExp.C

\*---------------------------------------------------------------------------*/

#ifndef regExp_H
#define regExp_H

#include <regex.h>
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class string;
template<class T> class List;

/*---------------------------------------------------------------------------*\
                          Class regExp Declaration
\*---------------------------------------------------------------------------*/

class regExp
{
    // Private data

        //- Precompiled regular expression
        mutable regex_t* preg_;

    // Private member functions

        //- Compile into a regular expression
        void compile(const char*) const;

        //- Disallow default bitwise copy construct
        regExp(const regExp&);

        //- Disallow default bitwise assignment
        void operator=(const regExp&);

public:

        //- Is character a regular expression meta-character?
        //  any character: '.' \n
        //  quantifiers: '*', '+', '?' \n
        //  grouping: '(', '|', ')' \n
        //  range: '[', ']' \n
        //
        //  Don't bother checking for '{digit}' bounds
        inline static bool meta(char c)
        {
            return
            (
                (c == '.')                           // any character
             || (c == '*' || c == '+' || c == '?')   // quantifiers
             || (c == '(' || c == ')' || c == '|')   // grouping/branching
             || (c == '[' || c == ']')               // range
            );
        }


    // Constructors

        //- Construct null
        regExp();

        //- Construct from character array
        regExp(const char*);

        //- Construct from std::string (or string)
        regExp(const std::string&);

    // Destructor

        ~regExp();


    // Member functions

        //- Is the precompiled expression set?
        inline bool exists() const
        {
            return preg_ ? true : false;
        }

        //- Return the number of (groups)
        inline int ngroups() const
        {
            return preg_ ? preg_->re_nsub : 0;
        }

        //- Release precompiled expression.
        //  Returns true if precompiled expression existed before clear
        bool clear() const;

        //- Return true if it matches, partial matches are optional
        bool match(const std::string&, bool partial=false) const;

        //- Return true if it matches and sets the sub-groups matched
        bool match(const string&, List<string>& groups) const;


    // Member Operators

        //- Assign from a string and compile regular expression
        void operator=(const std::string&);

        //- Assign from a character array and compile regular expression
        void operator=(const char*);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
