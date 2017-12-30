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
            if (iter.second)
            {
                // Emit warning, but only if version (non-zero) was provided
                std::cerr
                    << "--> FOAM IOWarning :" << nl
                    << "    Found [v" << iter.second << "] '"
                    << iter.first << "' instead of '"
                    << keyword.c_str() << "' in dictionary \""
                    << name().c_str() << "\" "
                    << nl
                    << std::endl;

                error::warnAboutAge("keyword", iter.second);
            }

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
