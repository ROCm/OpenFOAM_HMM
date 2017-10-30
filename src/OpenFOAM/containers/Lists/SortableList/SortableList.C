/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "ListOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::SortableList<T>::SortableList()
{}


template<class T>
inline Foam::SortableList<T>::SortableList(const label size)
:
    List<T>(size)
{}


template<class T>
inline Foam::SortableList<T>::SortableList(const label size, const T& val)
:
    List<T>(size, val)
{}


template<class T>
inline Foam::SortableList<T>::SortableList(const SortableList<T>& lst)
:
    List<T>(lst),
    indices_(lst.indices())
{}


template<class T>
inline Foam::SortableList<T>::SortableList(SortableList<T>&& lst)
:
    List<T>(std::move(lst)),
    indices_(std::move(lst.indices_))
{}


template<class T>
Foam::SortableList<T>::SortableList(const UList<T>& values)
:
    List<T>(values)
{
    sort();
}


template<class T>
Foam::SortableList<T>::SortableList(List<T>&& values)
:
    List<T>(std::move(values))
{
    sort();
}


template<class T>
Foam::SortableList<T>::SortableList(const Xfer<List<T>>& values)
:
    List<T>(values)
{
    sort();
}


template<class T>
template<class InputIterator>
inline Foam::SortableList<T>::SortableList
(
    InputIterator begIter,
    InputIterator endIter
)
:
    List<T>(begIter, endIter)
{
    sort();
}


template<class T>
Foam::SortableList<T>::SortableList(std::initializer_list<T> values)
:
    List<T>(values)
{
    sort();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::SortableList<T>::clear()
{
    List<T>::clear();
    indices_.clear();
}


template<class T>
Foam::List<T>& Foam::SortableList<T>::shrink()
{
    indices_.clear();
    return static_cast<List<T>&>(*this);
}


template<class T>
void Foam::SortableList<T>::sort()
{
    Foam::sortedOrder(*this, indices_);

    List<T> lst(*this, indices_); // Copy with indices for mapping
    List<T>::transfer(lst);
}


template<class T>
void Foam::SortableList<T>::reverseSort()
{
    Foam::sortedOrder(*this, indices_, typename UList<T>::greater(*this));

    List<T> lst(*this, indices_); // Copy with indices for mapping
    List<T>::transfer(lst);
}


template<class T>
void Foam::SortableList<T>::swap(SortableList<T>& lst)
{
    List<T>::swap(lst);
    indices_.swap(lst.indices_);
}


template<class T>
Foam::Xfer<Foam::List<T>> Foam::SortableList<T>::xfer()
{
    return xferMoveTo<List<T>>(*this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline void Foam::SortableList<T>::operator=(const T& val)
{
    indices_.clear();
    UList<T>::operator=(val);
}


template<class T>
inline void Foam::SortableList<T>::operator=(const UList<T>& lst)
{
    indices_.clear();
    List<T>::operator=(lst);
}


template<class T>
inline void Foam::SortableList<T>::operator=(const SortableList<T>& lst)
{
    List<T>::operator=(lst);
    indices_ = lst.indices();
}


template<class T>
inline void Foam::SortableList<T>::operator=(List<T>&& lst)
{
    indices_.clear();
    List<T>::operator=(std::move(lst));
}


template<class T>
inline void Foam::SortableList<T>::operator=(SortableList<T>&& lst)
{
    clear();
    this->swap(lst);
}


template<class T>
inline void Foam::SortableList<T>::operator=(std::initializer_list<T> lst)
{
    List<T>::operator=(lst);
    sort();
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class T>
inline void Foam::Swap(SortableList<T>& a, SortableList<T>& b)
{
    a.swap(b);
}


// ************************************************************************* //
