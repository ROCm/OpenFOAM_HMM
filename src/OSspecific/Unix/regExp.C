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

\*---------------------------------------------------------------------------*/

#include <sys/types.h>
#include "regExp.H"
#include "label.H"
#include "string.H"
#include "List.H"
#include "IOstreams.H"


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::regExp::compile(const char* pat) const
{
    clear();
    preg_ = new regex_t;

    if (regcomp(preg_, pat, REG_EXTENDED) != 0)
    {
        FatalErrorIn
        (
            "regExp::compile(const char*)"
        )   << "Failed to compile regular expression '" << pat << "'"
            << exit(FatalError);
    }
}


void Foam::regExp::clear() const
{
    if (preg_)
    {
        regfree(preg_);
        delete preg_;
        preg_ = 0;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regExp::regExp()
:
    preg_(0)
{}


Foam::regExp::regExp(const string& pat)
:
    preg_(0)
{
    compile(pat.c_str());
}


Foam::regExp::regExp(const char* pat)
:
    preg_(0)
{
    compile(pat);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regExp::~regExp()
{
    clear();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

int Foam::regExp::ngroups() const
{
    return preg_ ? preg_->re_nsub : 0;
}


bool Foam::regExp::match
(
    const string& str,
    bool partialMatch
) const
{
    if (preg_ && str.size())
    {
        size_t nmatch = 1;
        regmatch_t pmatch[1];

        // match and also verify that the entire string was matched
        if
        (
            regexec(preg_, str.c_str(), nmatch, pmatch, 0) == 0
         &&
            (
                partialMatch
             || (pmatch[0].rm_so == 0 && pmatch[0].rm_eo == label(str.size()))
            )
        )
        {
            return true;
        }
    }

    return false;
}


bool Foam::regExp::match
(
    const string& str,
    List<string>& groups,
    bool partialMatch
) const
{
    if (preg_ && str.size())
    {
        size_t nmatch = ngroups() + 1;
        regmatch_t pmatch[nmatch];

        // match and also verify that the entire string was matched
        if
        (
            regexec(preg_, str.c_str(), nmatch, pmatch, 0) == 0
         &&
            (
                partialMatch
             || (pmatch[0].rm_so == 0 && pmatch[0].rm_eo == label(str.size()))
            )
        )
        {
            groups.setSize(ngroups());
            label groupI = 0;

            for (size_t matchI = 1; matchI < nmatch; matchI++)
            {
                if (pmatch[matchI].rm_so != -1 && pmatch[matchI].rm_eo != -1)
                {
                    groups[groupI] = str.substr
                    (
                        pmatch[matchI].rm_so,
                        pmatch[matchI].rm_eo - pmatch[matchI].rm_so
                    );
                }
                else
                {
                    groups[groupI].clear();
                }
                groupI++;
            }

            return true;
        }
    }

    groups.clear();
    return false;
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


void Foam::regExp::operator=(const string& pat)
{
    compile(pat.c_str());
}


void Foam::regExp::operator=(const char* pat)
{
    compile(pat);
}

// ************************************************************************* //
