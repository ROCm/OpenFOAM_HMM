/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
inline Foam::IndirectList<T>::IndirectList
(
    const UList<T>& values,
    const labelUList& addr
)
:
    IndirectListAddressing<labelList>(addr),
    UIndirectList<T>
    (
        values,
        IndirectListAddressing<labelList>::addressing()
    )
{}


template<class T>
inline Foam::IndirectList<T>::IndirectList
(
    const UList<T>& values,
    labelList&& addr
)
:
    IndirectListAddressing<labelList>(std::move(addr)),
    UIndirectList<T>
    (
        values,
        IndirectListAddressing<labelList>::addressing()
    )
{}


template<class T>
inline Foam::IndirectList<T>::IndirectList
(
    const UList<T>& values,
    const Foam::zero
)
:
    IndirectListAddressing<labelList>(labelList()),  // zero-size addressing
    UIndirectList<T>
    (
        values,
        IndirectListAddressing<labelList>::addressing()
    )
{}


template<class T>
inline Foam::IndirectList<T>::IndirectList(const IndirectList<T>& list)
:
    // Copy addressing
    IndirectListAddressing<labelList>(list.addressing()),
    UIndirectList<T>
    (
        list.values(),
        IndirectListAddressing<labelList>::addressing()
    )
{}


template<class T>
inline Foam::IndirectList<T>::IndirectList(IndirectList<T>&& list)
:
    // Move addressing
    IndirectListAddressing<labelList>(std::move(list.addressing())),
    UIndirectList<T>
    (
        list.values(),
        IndirectListAddressing<labelList>::addressing()
    )
{}


template<class T>
inline Foam::IndirectList<T>::IndirectList(const UIndirectList<T>& list)
:
    // Copy addressing
    IndirectListAddressing<labelList>(list.addressing()),
    UIndirectList<T>
    (
        list.values(),
        IndirectListAddressing<labelList>::addressing()
    )
{}


// ************************************************************************* //
