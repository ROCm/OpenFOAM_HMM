/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "dictionary.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

static void warnAboutAge(const int oldVersion)
{
    if (oldVersion < 1000)
    {
        // Emit warning
        std::cerr
            << "    This keyword is considered to be VERY old!\n"
            << std::endl;
    }
#if (OPENFOAM_PLUS > 1600)
    else if (OPENFOAM_PLUS > oldVersion)
    {
        const int months =
        (
            // YYMM -> months
            (12 * (OPENFOAM_PLUS/100) + (OPENFOAM_PLUS % 100))
          - (12 * (oldVersion/100) + (oldVersion % 100))
        );

        std::cerr
            << "    This keyword is deemed to be " << months
            << " months old.\n"
            << std::endl;
    }
#endif
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dictionary::const_searcher Foam::dictionary::csearchCompat
(
    const word& keyword,
    std::initializer_list<std::pair<const char*,int>> compat,
    bool recursive,
    bool patternMatch
) const
{
    const_searcher finder(csearch(keyword, recursive,patternMatch));

    if (finder.found())
    {
        return finder;
    }

    for (const std::pair<const char*,int>& iter : compat)
    {
        finder = csearch(word::validate(iter.first), recursive,patternMatch);

        if (finder.found())
        {
            // Emit warning
            std::cerr
                << "--> FOAM IOWarning :" << nl
                << "    Found [v" << iter.second << "] '"
                << iter.first << "' instead of '"
                << keyword.c_str() << "' in dictionary \""
                << name().c_str() << "\" "
                << nl
                << std::endl;

            warnAboutAge(iter.second);

            break;
        }
    }

    return finder;
}


bool Foam::dictionary::foundCompat
(
    const word& keyword,
    std::initializer_list<std::pair<const char*,int>> compat,
    bool recursive,
    bool patternMatch
) const
{
    return csearchCompat(keyword, compat, recursive,patternMatch).found();
}


const Foam::entry* Foam::dictionary::lookupEntryPtrCompat
(
    const word& keyword,
    std::initializer_list<std::pair<const char*,int>> compat,
    bool recursive,
    bool patternMatch
) const
{
    return csearchCompat(keyword, compat, recursive,patternMatch).ptr();
}


const Foam::entry& Foam::dictionary::lookupEntryCompat
(
    const word& keyword,
    std::initializer_list<std::pair<const char*,int>> compat,
    bool recursive,
    bool patternMatch
) const
{
    const const_searcher
        finder(csearchCompat(keyword, compat, recursive,patternMatch));

    if (!finder.found())
    {
        FatalIOErrorInFunction
        (
            *this
        )   << "keyword " << keyword << " is undefined in dictionary "
            << name()
            << exit(FatalIOError);
    }

    return finder.ref();
}


Foam::ITstream& Foam::dictionary::lookupCompat
(
    const word& keyword,
    std::initializer_list<std::pair<const char*,int>> compat,
    bool recursive,
    bool patternMatch
) const
{
    return lookupEntryCompat(keyword, compat, recursive,patternMatch).stream();
}


// ************************************************************************* //
