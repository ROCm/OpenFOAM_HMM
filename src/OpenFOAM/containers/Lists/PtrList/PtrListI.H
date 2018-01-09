/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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

#include "autoPtr.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
inline Foam::PtrList<T>::PtrList(const label nElem)
:
    UPtrList<T>(nElem)
{}


template<class T>
inline Foam::PtrList<T>::PtrList(PtrList<T>&& lst)
:
    UPtrList<T>(std::move(lst))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline void Foam::PtrList<T>::resize(const label newLen)
{
    this->setSize(newLen);
}


template<class T>
inline void Foam::PtrList<T>::append(const autoPtr<T>& aptr)
{
    return UPtrList<T>::append(const_cast<autoPtr<T>&>(aptr).ptr());
}


template<class T>
inline void Foam::PtrList<T>::append(const tmp<T>& tptr)
{
    return UPtrList<T>::append(const_cast<tmp<T>&>(tptr).ptr());
}


template<class T>
inline bool Foam::PtrList<T>::set(const label i) const
{
    return this->ptrs_[i] != nullptr;
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrList<T>::set(const label i, T* ptr)
{
    autoPtr<T> old(this->ptrs_[i]);
    this->ptrs_[i] = ptr;
    return old;
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrList<T>::set
(
    const label i,
    const autoPtr<T>& aptr
)
{
    return set(i, const_cast<autoPtr<T>&>(aptr).ptr());
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrList<T>::set
(
    const label i,
    const tmp<T>& tptr
)
{
    return set(i, const_cast<tmp<T>&>(tptr).ptr());
}


template<class T>
inline void Foam::PtrList<T>::transfer(PtrList<T>& lst)
{
    this->clear();
    UPtrList<T>::swap(lst);
}


template<class T>
inline Foam::Xfer<Foam::PtrList<T>> Foam::PtrList<T>::xfer()
{
    return xferMove(*this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline void Foam::PtrList<T>::operator=(PtrList<T>&& lst)
{
    this->clear();
    this->swap(lst);
}


// ************************************************************************* //
