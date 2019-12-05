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

#include "dictionary.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Should issue warning if there is +ve versioning (+ve version number)
// and the this version number is older than the current OpenFOAM version
// as conveyed by the OPENFOAM compiler define.

static inline constexpr bool shouldWarnVersion(const int version)
{
    return (version > 0 && version < OPENFOAM);
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dictionary::const_searcher Foam::dictionary::csearchCompat
(
    const word& keyword,
    std::initializer_list<std::pair<const char*,int>> compat,
    enum keyType::option matchOpt
) const
{
    const_searcher finder(csearch(keyword, matchOpt));

    if (finder.good())
    {
        return finder;
    }

    for (const std::pair<const char*,int>& iter : compat)
    {
        finder = csearch(word::validate(iter.first), matchOpt);

        if (finder.good())
        {
            // Want single warning (on master), but guard with parRun to avoid
            // Pstream::master() when Pstream has not yet been initialized
            if
            (
                shouldWarnVersion(iter.second)
             && (Pstream::parRun() ? Pstream::master() : true)
            )
            {
                std::cerr
                    << "--> FOAM IOWarning :" << nl
                    << "    Found [v" << iter.second << "] '"
                    << iter.first << "' entry instead of '"
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
    enum keyType::option matchOpt
) const
{
    return csearchCompat(keyword, compat, matchOpt).good();
}


const Foam::entry* Foam::dictionary::findCompat
(
    const word& keyword,
    std::initializer_list<std::pair<const char*,int>> compat,
    enum keyType::option matchOpt
) const
{
    return csearchCompat(keyword, compat, matchOpt).ptr();
}


const Foam::entry& Foam::dictionary::lookupEntryCompat
(
    const word& keyword,
    std::initializer_list<std::pair<const char*,int>> compat,
    enum keyType::option matchOpt
) const
{
    const const_searcher finder(csearchCompat(keyword, compat, matchOpt));

    if (!finder.good())
    {
        FatalIOErrorInFunction(*this)
            << "Entry '" << keyword << "' not found in dictionary "
            << name()
            << exit(FatalIOError);
    }

    return finder.ref();
}


Foam::ITstream& Foam::dictionary::lookupCompat
(
    const word& keyword,
    std::initializer_list<std::pair<const char*,int>> compat,
    enum keyType::option matchOpt
) const
{
    return lookupEntryCompat(keyword, compat, matchOpt).stream();
}


// ************************************************************************* //
