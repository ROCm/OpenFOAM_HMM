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

#include "Enum.H"
#include "dictionary.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EnumType>
Foam::label Foam::Enum<EnumType>::getIndex(const word& enumName) const
{
    const label n = size();
    for (label idx=0; idx < n; ++idx)
    {
        if (names_[idx] == enumName)
        {
            return idx;
        }
    }

    return -1;
}


template<class EnumType>
Foam::label Foam::Enum<EnumType>::getIndex(const EnumType e) const
{
    const int val = int(e);

    const label n = size();
    for (label idx=0; idx < n; ++idx)
    {
        if (values_[idx] == val)
        {
            return idx;
        }
    }

    return -1;
}


template<class EnumType>
EnumType Foam::Enum<EnumType>::getEnum(const word& enumName) const
{
    const label idx = getIndex(enumName);

    if (idx < 0)
    {
        FatalErrorInFunction
            << enumName << " is not in enumeration: "
            << names_ << exit(FatalError);
    }

    return EnumType(values_[idx]);
}


template<class EnumType>
const Foam::word& Foam::Enum<EnumType>::getName(const EnumType e) const
{
    const label idx = getIndex(e);

    if (idx < 0)
    {
        return word::null;
    }

    return names_[idx];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EnumType>
Foam::Enum<EnumType>::Enum
(
    std::initializer_list<std::pair<EnumType, word>> lst
)
:
    names_(lst.size()),
    values_(lst.size())
{
    int idx = 0;
    for (const auto& pair : lst)
    {
        names_[idx]  = pair.second;
        values_[idx] = int(pair.first);

        ++idx;
    }
}


template<class EnumType>
Foam::Enum<EnumType>::Enum
(
    const EnumType start,
    std::initializer_list<word> lst
)
:
    names_(lst.size()),
    values_(lst.size())
{
    int val = int(start);

    int idx = 0;
    for (const auto& key : lst)
    {
        names_[idx]  = key;
        values_[idx] = val;

        ++val;
        ++idx;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EnumType>
Foam::List<Foam::word> Foam::Enum<EnumType>::sortedToc() const
{
    wordList lst(names_);
    Foam::sort(lst);

    return lst;
}


template<class EnumType>
EnumType Foam::Enum<EnumType>::read(Istream& is) const
{
    const word enumName(is);
    const label idx = getIndex(enumName);

    if (idx < 0)
    {
        FatalIOErrorInFunction(is)
            << enumName << " is not in enumeration: "
            << names_ << nl
            << exit(FatalIOError);
    }

    return EnumType(values_[idx]);
}


template<class EnumType>
void Foam::Enum<EnumType>::write(const EnumType e, Ostream& os) const
{
    const label idx = getIndex(e);
    if (idx >= 0)
    {
        os << names_[idx];
    }
}


template<class EnumType>
EnumType Foam::Enum<EnumType>::lookup
(
    const word& key,
    const dictionary& dict
) const
{
    const word enumName(dict.lookup(key));
    const label idx = getIndex(enumName);

    if (idx < 0)
    {
        FatalIOErrorInFunction(dict)
            << enumName << " is not in enumeration: "
            << names_ << nl
            << exit(FatalIOError);
    }

    return EnumType(values_[idx]);
}


template<class EnumType>
EnumType Foam::Enum<EnumType>::lookupOrDefault
(
    const word& key,
    const dictionary& dict,
    const EnumType deflt
) const
{
    if (dict.found(key))
    {
        return lookup(key, dict);
    }

    return deflt;
}


template<class EnumType>
EnumType Foam::Enum<EnumType>::lookupOrFailsafe
(
    const word& key,
    const dictionary& dict,
    const EnumType deflt
) const
{
    if (dict.found(key))
    {
        const word enumName(dict.lookup(key));
        const label idx = getIndex(enumName);

        if (idx >= 0)
        {
            return EnumType(values_[idx]);
        }

        IOWarningInFunction(dict)
            << "bad " << key <<" specifier " << enumName
            << " using " << getName(deflt) << endl;
    }

    return deflt;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class EnumType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const Enum<EnumType>& wrapped
)
{
    return wrapped.names().writeList(os, 10);
}


// ************************************************************************* //
