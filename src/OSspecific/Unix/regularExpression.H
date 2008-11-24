/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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
    Foam::regularExpression

Description
    Wrapper around regular expressions.

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef regularExpression_H
#define regularExpression_H

#include <sys/types.h>
#include <regex.h>
#include "string.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class regularExpression Declaration
\*---------------------------------------------------------------------------*/

class regularExpression
{
    // Private data

        //- Precompiled regular expression
        regex_t* preg_;


    // Private member functions

        //- Disallow default bitwise copy construct
        regularExpression(const regularExpression&);

        //- Disallow default bitwise assignment
        void operator=(const regularExpression&);

public:


    // Constructors

        //- Construct from string
        inline regularExpression(const string& s)
        {
            preg_ = new regex_t;

            if (regcomp(preg_, s.c_str(), REG_EXTENDED) != 0)
            {
                FatalErrorIn
                (
                    "regularExpression::regularExpression(const char*)"
                )   << "Failed to compile regular expression " << s
                    << exit(FatalError);
            }
        }


    // Destructor

        //- Construct from string
        inline ~regularExpression()
        {
            if (preg_)
            {
                regfree(preg_);
                delete preg_;
            }
        }


    // Member functions

        //- Matches?
        inline bool matches(const string& s) const
        {
            size_t nmatch = 1;
            regmatch_t pmatch[1];

            int errVal = regexec(preg_, s.c_str(), nmatch, pmatch, 0);

            return (errVal == 0 && pmatch[0].rm_eo == label(s.size()));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
