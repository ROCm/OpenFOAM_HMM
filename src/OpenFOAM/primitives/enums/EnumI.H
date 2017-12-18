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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EnumType>
inline Foam::label Foam::Enum<EnumType>::size() const
{
    return names_.size();
}


template<class EnumType>
inline const Foam::wordList& Foam::Enum<EnumType>::names() const
{
    return names_;
}


template<class EnumType>
inline const Foam::wordList& Foam::Enum<EnumType>::toc() const
{
    return names_;
}


template<class EnumType>
inline const Foam::List<int>& Foam::Enum<EnumType>::values() const
{
    return values_;
}


template<class EnumType>
inline bool Foam::Enum<EnumType>::found(const word& enumName) const
{
    return getIndex(enumName) >= 0;
}


template<class EnumType>
inline bool Foam::Enum<EnumType>::found(const EnumType e) const
{
    return getIndex(e) >= 0;
}


template<class EnumType>
inline bool Foam::Enum<EnumType>::hasEnum(const word& enumName) const
{
    return getIndex(enumName) >= 0;
}


template<class EnumType>
inline bool Foam::Enum<EnumType>::hasName(const EnumType e) const
{
    return getIndex(e) >= 0;
}


template<class EnumType>
inline Foam::Ostream& Foam::Enum<EnumType>::writeList
(
    Ostream& os,
    const label shortListLen
) const
{
    return names_.writeList(os, shortListLen);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EnumType>
inline EnumType Foam::Enum<EnumType>::operator[]
(
    const word& enumName
) const
{
    return getEnum(enumName);
}


template<class EnumType>
inline const Foam::word& Foam::Enum<EnumType>::operator[]
(
    const EnumType e
) const
{
    return getName(e);
}


template<class EnumType>
inline EnumType Foam::Enum<EnumType>::operator()
(
    const word& enumName,
    const EnumType deflt
) const
{
    const label idx = getIndex(enumName);

    if (idx >= 0)
    {
        return EnumType(values_[idx]);
    }

    return deflt;
}


// ************************************************************************* //
