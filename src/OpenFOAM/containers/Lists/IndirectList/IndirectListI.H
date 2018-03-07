/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::IndirectListAddressing::IndirectListAddressing
(
    const labelUList& addr
)
:
    addressing_(addr)
{}


inline Foam::IndirectListAddressing::IndirectListAddressing
(
    List<label>&& addr
)
:
    addressing_(std::move(addr))
{}


template<class T>
inline Foam::IndirectList<T>::IndirectList
(
    const UList<T>& values,
    const labelUList& addr
)
:
    IndirectListAddressing(addr),
    UIndirectList<T>
    (
        values,
        IndirectListAddressing::addressing()
    )
{}


template<class T>
inline Foam::IndirectList<T>::IndirectList
(
    const UList<T>& values,
    List<label>&& addr
)
:
    IndirectListAddressing(std::move(addr)),
    UIndirectList<T>
    (
        values,
        IndirectListAddressing::addressing()
    )
{}


template<class T>
inline Foam::IndirectList<T>::IndirectList
(
    const IndirectList<T>& list
)
:
    IndirectListAddressing(list.addressing()),  // Copy addressing
    UIndirectList<T>
    (
        list.completeList(),
        IndirectListAddressing::addressing()
    )
{}


template<class T>
inline Foam::IndirectList<T>::IndirectList
(
    const UIndirectList<T>& list
)
:
    IndirectListAddressing(list.addressing()),  // Copy addressing
    UIndirectList<T>
    (
        list.completeList(),
        IndirectListAddressing::addressing()
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::List<Foam::label>&
Foam::IndirectListAddressing::addressing() const
{
    return addressing_;
}


inline void Foam::IndirectListAddressing::resetAddressing
(
    const labelUList& addr
)
{
    addressing_ = addr;
}


inline void Foam::IndirectListAddressing::resetAddressing
(
    List<label>&& addr
)
{
    addressing_.transfer(addr);
}


// ************************************************************************* //
