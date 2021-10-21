/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::SubList<T>::SubList
(
    const UList<T>& list
)
:
    UList<T>(const_cast<T*>(list.cdata()), list.size())
{}


template<class T>
template<unsigned N>
inline Foam::SubList<T>::SubList
(
    const FixedList<T, N>& list
)
:
    UList<T>(const_cast<T*>(list.cdata()), static_cast<label>(N))
{}


template<class T>
inline Foam::SubList<T>::SubList
(
    const UList<T>& list,
    const label subSize
)
:
    UList<T>(const_cast<T*>(list.cdata()), subSize)
{
    #ifdef FULLDEBUG
    list.checkSize(subSize);
    #endif
}


template<class T>
inline Foam::SubList<T>::SubList
(
    const UList<T>& list,
    const label subSize,
    const label startIndex
)
:
    UList<T>(const_cast<T*>(list.cdata() + startIndex), subSize)
{
    #ifdef FULLDEBUG
    list.checkRange(startIndex, subSize);
    #endif
}


template<class T>
inline Foam::SubList<T>::SubList
(
    const UList<T>& list,
    const labelRange& range
)
:
    SubList<T>(range.subset0(list.size()), list)
{
    #ifdef FULLDEBUG
    // The subset0() above always produces valid ranges but want to check
    // that the input itself was valid
    list.checkRange(range.start(), range.size());
    #endif
}


template<class T>
inline Foam::SubList<T>::SubList
(
    const labelRange& range,
    const UList<T>& list
)
:
    SubList<T>(list, range.size(), range.start())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline const Foam::SubList<T>& Foam::SubList<T>::null()
{
    return NullObjectRef<SubList<T>>();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline Foam::SubList<T>::operator const Foam::List<T>&() const
{
    return *reinterpret_cast<const List<T>*>(this);
}


template<class T>
inline void Foam::SubList<T>::operator=(const SubList<T>& list)
{
    UList<T>::deepCopy(list);
}


template<class T>
inline void Foam::SubList<T>::operator=(const UList<T>& list)
{
    UList<T>::deepCopy(list);
}


template<class T>
template<class Addr>
inline void Foam::SubList<T>::operator=(const IndirectListBase<T, Addr>& list)
{
    UList<T>::deepCopy(list);
}


template<class T>
inline void Foam::SubList<T>::operator=(const T& val)
{
    UList<T>::operator=(val);
}


template<class T>
inline void Foam::SubList<T>::operator=(const Foam::zero)
{
    UList<T>::operator=(Zero);
}


// ************************************************************************* //
