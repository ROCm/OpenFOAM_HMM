/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

#include "regExpCxx.H"
#include "debug.H"
#include "error.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::regExpCxx::grammar(Foam::debug::optimisationSwitch("regExpCxx", 0));


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

static std::string error_string(const std::regex_error& err)
{
    switch (err.code())
    {
        case std::regex_constants::error_collate :
            return "invalid collating element name";
            break;

        case std::regex_constants::error_ctype :
            return "invalid character class name";
            break;

        case std::regex_constants::error_escape :
            return "invalid escaped character or a trailing escape";
            break;

        case std::regex_constants::error_backref :
            return "invalid back reference";
            break;

        case std::regex_constants::error_brack :
            return "mismatched [ and ]";
            break;

        case std::regex_constants::error_paren :
            return "mismatched ( and )";
            break;

        case std::regex_constants::error_brace :
            return "mismatched { and }";
            break;

        case std::regex_constants::error_badbrace :
            return "invalid range in a {..}";
            break;

        case std::regex_constants::error_range :
            return "invalid [..] character range";
            break;

        case std::regex_constants::error_space :
            return "memory error";
            break;

        case std::regex_constants::error_badrepeat :
            return "bad '*?+{' repeat";
            break;

        case std::regex_constants::error_complexity :
            return "expression too complex";
            break;

        case std::regex_constants::error_stack :
            return "memory stack error";
            break;

        default:
            break;

    }

    return "";
}

} // End anonymous namespace


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::regExpCxx::set(const char* pattern, bool ignoreCase)
{
    clear();  // Also sets ok_ = false

    size_t len = (pattern ? strlen(pattern) : 0);

    // Avoid nullptr and zero-length patterns
    if (!len)
    {
        return false;
    }

    std::regex::flag_type flags = syntax();
    if (ignoreCase)
    {
        flags |= std::regex::icase;
    }

    const char* pat = pattern;

    // Has embedded ignore-case prefix?
    if (len >= 4 && !strncmp(pattern, "(?i)", 4))
    {
        flags |= std::regex::icase;
        pat += 4;
        len -= 4;
    }

    if (len)
    {
        try
        {
            re_.assign(pat, flags);
            ok_ = true;
        }
        catch (const std::regex_error& err)
        {
            FatalErrorInFunction
                << "Failed to compile regular expression '"
                << pattern << "'" << nl
                << err.what() << ": " << error_string(err).c_str() << nl
                << exit(FatalError);
        }
    }

    return ok_;
}


bool Foam::regExpCxx::set(const std::string& pattern, bool ignoreCase)
{
    clear();  // Also sets ok_ = false

    auto len = pattern.size();

    // Avoid zero-length patterns
    if (!len)
    {
        return false;
    }

    std::regex::flag_type flags = syntax();
    if (ignoreCase)
    {
        flags |= std::regex::icase;
    }

    auto pat = pattern.begin();

    // Has embedded ignore-case prefix?
    if (len >= 4 && !pattern.compare(0, 4, "(?i)"))
    {
        flags |= std::regex::icase;
        pat += 4;
        len -= 4;
    }

    if (len)
    {
        try
        {
            re_.assign(pat, pattern.end(), flags);
            ok_ = true;
        }
        catch (const std::regex_error& err)
        {
            FatalErrorInFunction
                << "Failed to compile regular expression '"
                << pattern.c_str() << "'" << nl
                << err.what() << ": " << error_string(err).c_str() << nl
                << exit(FatalError);
        }
    }

    return ok_;
}


// ************************************************************************* //
