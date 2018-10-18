/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenCFD Ltd.
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
    return keys_.size();
}


template<class EnumType>
inline const Foam::wordList& Foam::Enum<EnumType>::names() const
{
    return keys_;
}


template<class EnumType>
inline const Foam::wordList& Foam::Enum<EnumType>::toc() const
{
    return keys_;
}


template<class EnumType>
inline const Foam::List<int>& Foam::Enum<EnumType>::values() const
{
    return vals_;
}


template<class EnumType>
inline Foam::label Foam::Enum<EnumType>::find(const word& enumName) const
{
    return keys_.find(enumName);
}


template<class EnumType>
inline Foam::label Foam::Enum<EnumType>::find(const EnumType e) const
{
    return vals_.find(int(e));
}


template<class EnumType>
inline bool Foam::Enum<EnumType>::found(const word& enumName) const
{
    return keys_.found(enumName);
}


template<class EnumType>
inline bool Foam::Enum<EnumType>::found(const EnumType e) const
{
    return vals_.found(int(e));
}


template<class EnumType>
inline const Foam::word& Foam::Enum<EnumType>::get(const EnumType e) const
{
    const label idx = find(e);

    if (idx < 0)
    {
        return word::null;
    }

    return keys_[idx];
}


template<class EnumType>
inline Foam::Ostream& Foam::Enum<EnumType>::writeList
(
    Ostream& os,
    const label shortLen
) const
{
    return keys_.writeList(os, shortLen);
}


template<class EnumType>
inline void Foam::Enum<EnumType>::write(const EnumType e, Ostream& os) const
{
    const label idx = find(e);

    if (idx >= 0)
    {
        os << keys_[idx];
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EnumType>
inline EnumType Foam::Enum<EnumType>::operator[]
(
    const word& enumName
) const
{
    return get(enumName);
}


template<class EnumType>
inline const Foam::word& Foam::Enum<EnumType>::operator[]
(
    const EnumType e
) const
{
    return get(e);
}


template<class EnumType>
inline EnumType Foam::Enum<EnumType>::operator()
(
    const word& enumName,
    const EnumType defaultValue
) const
{
    const label idx = find(enumName);

    if (idx >= 0)
    {
        return EnumType(vals_[idx]);
    }

    return defaultValue;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class EnumType>
inline Foam::Ostream& Foam::operator<<(Ostream& os, const Enum<EnumType>& list)
{
    return list.names().writeList(os, 0);
}


// ************************************************************************* //
