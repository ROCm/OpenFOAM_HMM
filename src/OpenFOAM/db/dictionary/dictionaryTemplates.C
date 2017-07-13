/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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
#include "primitiveEntry.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Compare>
Foam::wordList Foam::dictionary::sortedToc(const Compare& comp) const
{
    return hashedEntries_.sortedToc(comp);
}


template<class T>
T Foam::dictionary::lookupType
(
    const word& keyword,
    bool recursive,
    bool patternMatch
) const
{
    auto finder = csearch(keyword, recursive, patternMatch);

    if (!finder.found())
    {
        FatalIOErrorInFunction
        (
            *this
        )   << "keyword " << keyword << " is undefined in dictionary "
            << name()
            << exit(FatalIOError);
    }

    return pTraits<T>(finder.ptr()->stream());
}


template<class T>
T Foam::dictionary::lookupOrDefault
(
    const word& keyword,
    const T& deflt,
    bool recursive,
    bool patternMatch
) const
{
    auto finder = csearch(keyword, recursive, patternMatch);

    if (finder.found())
    {
        return pTraits<T>(finder.ptr()->stream());
    }

    if (writeOptionalEntries)
    {
        IOInfoInFunction(*this)
            << "Optional entry '" << keyword << "' is not present,"
            << " returning the default value '" << deflt << "'"
            << endl;
    }

    return deflt;
}


template<class T>
T Foam::dictionary::lookupOrAddDefault
(
    const word& keyword,
    const T& deflt,
    bool recursive,
    bool patternMatch
)
{
    auto finder = csearch(keyword, recursive, patternMatch);

    if (finder.found())
    {
        return pTraits<T>(finder.ptr()->stream());
    }

    if (writeOptionalEntries)
    {
        IOInfoInFunction(*this)
            << "Optional entry '" << keyword << "' is not present,"
            << " adding and returning the default value '" << deflt << "'"
            << endl;
    }

    add(new primitiveEntry(keyword, deflt));
    return deflt;
}


template<class T>
bool Foam::dictionary::readIfPresent
(
    const word& keyword,
    T& val,
    bool recursive,
    bool patternMatch
) const
{
    auto finder = csearch(keyword, recursive, patternMatch);

    if (finder.found())
    {
        finder.ptr()->stream() >> val;
        return true;
    }

    if (writeOptionalEntries)
    {
        IOInfoInFunction(*this)
            << "Optional entry '" << keyword << "' is not present,"
            << " the default value '" << val << "' will be used."
            << endl;
    }

    return false;
}


template<class T>
void Foam::dictionary::add(const keyType& k, const T& v, bool overwrite)
{
    add(new primitiveEntry(k, v), overwrite);
}


template<class T>
void Foam::dictionary::set(const keyType& k, const T& v)
{
    set(new primitiveEntry(k, v));
}


// ************************************************************************* //
