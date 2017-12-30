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

// * * * * * * * * * * * * * * * helper methods  * * * * * * * * * * * * * * //

template<class IteratorType, class TableType>
inline IteratorType Foam::HashTableCore::iterator_begin
(
    TableType& table
)
{
    return IteratorType(table.begin());
}


template<class IteratorType, class TableType>
inline IteratorType Foam::HashTableCore::iterator_cbegin
(
    const TableType& table
)
{
    return IteratorType(table.cbegin());
}


template<class IteratorType>
inline const IteratorType& Foam::HashTableCore::iterator_end()
{
    return *reinterpret_cast<const IteratorType*>(nullObjectPtr);
}


template<class IteratorType>
inline const IteratorType& Foam::HashTableCore::iterator_cend()
{
    return *reinterpret_cast<const IteratorType*>(nullObjectPtr);
}


// * * * * * * * * * * * * * const iterator pair * * * * * * * * * * * * * * //

template<class IteratorType, class TableType>
inline Foam::HashTableCore::const_iterator_pair<IteratorType, TableType>
::const_iterator_pair
(
    const TableType& tbl
)
:
    size_(tbl.size()),
    iter_(tbl.begin())
{}


template<class IteratorType, class TableType>
inline Foam::label
Foam::HashTableCore::const_iterator_pair<IteratorType, TableType>::size() const
{
    return size_;
}


template<class IteratorType, class TableType>
inline bool
Foam::HashTableCore::const_iterator_pair<IteratorType, TableType>::empty() const
{
    return !size_;
}


template<class IteratorType, class TableType>
inline IteratorType Foam::HashTableCore::const_iterator_pair
<
    IteratorType,
    TableType
>::begin() const
{
    return iter_;
}


template<class IteratorType, class TableType>
inline IteratorType Foam::HashTableCore::const_iterator_pair
<
    IteratorType,
    TableType
>::cbegin() const
{
    return iter_;
}


template<class IteratorType, class TableType>
inline const IteratorType& Foam::HashTableCore::const_iterator_pair
<
    IteratorType,
    TableType
>::end() const
{
    return HashTableCore::iterator_cend<IteratorType>();
}


template<class IteratorType, class TableType>
inline const IteratorType& Foam::HashTableCore::const_iterator_pair
<
    IteratorType,
    TableType
>::cend() const
{
    return HashTableCore::iterator_cend<IteratorType>();
}


// ************************************************************************* //
