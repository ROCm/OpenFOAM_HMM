/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "regExpPosix.H"
#include "SubStrings.H"
#include "error.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::regExpPosix::grammar(0);


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Matched entire length
static inline bool fullMatch(const regmatch_t& m, const regoff_t len)
{
    return (m.rm_so == 0 && m.rm_eo == len);
}

} // End anonymous namespace


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::regExpPosix::clear()
{
    if (preg_)
    {
        regfree(preg_);
        delete preg_;
        preg_ = nullptr;

        return true;
    }

    return false;
}


bool Foam::regExpPosix::set(const char* pattern, bool ignoreCase)
{
    clear();

    // Avoid nullptr and zero-length patterns
    if (pattern && *pattern)
    {
        int cflags = REG_EXTENDED;
        if (ignoreCase)
        {
            cflags |= REG_ICASE;
        }

        const char* pat = pattern;

        // Check for embedded prefix for ignore-case
        // this is the only embedded prefix we support
        // - a simple check is sufficient
        if (!strncmp(pattern, "(?i)", 4))
        {
            cflags |= REG_ICASE;
            pat += 4;

            // avoid zero-length patterns
            if (!*pat)
            {
                return false;
            }
        }

        preg_ = new regex_t;
        int err = regcomp(preg_, pat, cflags);

        if (err != 0)
        {
            char errbuf[200];
            regerror(err, preg_, errbuf, sizeof(errbuf));

            FatalErrorInFunction
                << "Failed to compile regular expression '" << pattern << "'"
                << nl << errbuf
                << exit(FatalError);
        }

        return true;
    }

    return false;  // Was cleared and nothing was set
}


bool Foam::regExpPosix::set(const std::string& pattern, bool ignoreCase)
{
    return set(pattern.c_str(), ignoreCase);
}


std::string::size_type Foam::regExpPosix::find(const std::string& text) const
{
    if (preg_ && !text.empty())
    {
        const size_t nmatch = 1;
        regmatch_t pmatch[1];

        if (regexec(preg_, text.c_str(), nmatch, pmatch, 0) == 0)
        {
            return pmatch[0].rm_so;
        }
    }

    return std::string::npos;
}


bool Foam::regExpPosix::match(const std::string& text) const
{
    const auto len = text.size();

    if (preg_ && len)
    {
        const size_t nmatch = 1;
        regmatch_t pmatch[1];

        // Verify that the entire string was matched
        // - [0] is the entire match result
        return
        (
            regexec(preg_, text.c_str(), nmatch, pmatch, 0) == 0
         && fullMatch(pmatch[0], len)
        );
    }

    return false;
}


bool Foam::regExpPosix::match
(
    const std::string& text,
    SubStrings<std::string>& matches
) const
{
    matches.clear();

    const auto len = text.size();
    if (preg_ && len)
    {
        const size_t nmatch = ngroups() + 1;
        regmatch_t pmatch[nmatch];

        // Verify that the entire string was matched
        // - [0] is the entire match result
        // - [1..] are the match groups (1..)
        if
        (
            regexec(preg_, text.c_str(), nmatch, pmatch, 0) != 0
         || !fullMatch(pmatch[0], len)
        )
        {
            return false;
        }

        matches.reserve(nmatch);

        for (size_t matchi = 0; matchi < nmatch; ++matchi)
        {
            const auto& mat = pmatch[matchi];

            if (mat.rm_so != -1 && mat.rm_eo != -1)
            {
                matches.append
                (
                    text.cbegin() + mat.rm_so,
                    text.cbegin() + mat.rm_eo
                );
            }
            else
            {
                // This may be misleading...
                matches.append(text.cbegin(), text.cbegin());
            }
        }

        return true;
    }

    return false;
}


// ************************************************************************* //
