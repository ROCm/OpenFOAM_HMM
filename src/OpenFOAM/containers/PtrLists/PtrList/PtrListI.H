/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2019 OpenCFD Ltd.
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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class T>
inline void Foam::PtrList<T>::free()
{
    (this->ptrs_).free(); // free old pointers
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
inline constexpr Foam::PtrList<T>::PtrList() noexcept
:
    UPtrList<T>()
{}


template<class T>
inline Foam::PtrList<T>::PtrList(const label len)
:
    UPtrList<T>(len)
{}


template<class T>
inline Foam::PtrList<T>::PtrList(const PtrList<T>& list)
:
    UPtrList<T>(list.ptrs_.clone())
{}


template<class T>
inline Foam::PtrList<T>::PtrList(PtrList<T>&& list)
:
    UPtrList<T>(std::move(list))
{}


template<class T>
inline Foam::PtrList<T>::PtrList(UList<T*>& list)
:
    UPtrList<T>(list.size())
{
    // Take ownership of the pointer
    forAll(list, i)
    {
        set(i, list[i]);
        list[i] = nullptr;
    }
}


template<class T>
template<class CloneArg>
inline Foam::PtrList<T>::PtrList
(
    const PtrList<T>& list,
    const CloneArg& cloneArg
)
:
    UPtrList<T>(list.clone(cloneArg)())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline void Foam::PtrList<T>::clear()
{
    (this->ptrs_).free(); // free old pointers
    UPtrList<T>::clear();
}


template<class T>
inline void Foam::PtrList<T>::setSize(const label newLen)
{
    this->resize(newLen);
}


template<class T>
inline void Foam::PtrList<T>::append(T* ptr)
{
    UPtrList<T>::append(ptr);
}


template<class T>
inline void Foam::PtrList<T>::append(autoPtr<T>& aptr)
{
    UPtrList<T>::append(aptr.release());
}


template<class T>
inline void Foam::PtrList<T>::append(autoPtr<T>&& aptr)
{
    UPtrList<T>::append(aptr.release());
}


template<class T>
inline void Foam::PtrList<T>::append(const tmp<T>& tptr)
{
    UPtrList<T>::append(tptr.ptr());
}


template<class T>
inline const T* Foam::PtrList<T>::set(const label i) const
{
    return UPtrList<T>::set(i);
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrList<T>::set(const label i, T* ptr)
{
    return autoPtr<T>(UPtrList<T>::set(i, ptr));
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrList<T>::set(const label i, autoPtr<T>& aptr)
{
    return set(i, aptr.release());
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrList<T>::set(const label i, autoPtr<T>&& aptr)
{
    return set(i, aptr.release());
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrList<T>::set(const label i, const tmp<T>& tptr)
{
    return set(i, tptr.ptr());
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrList<T>::release(const label i)
{
    if (i < 0 || i >= this->size())
    {
        return nullptr;
    }

    return autoPtr<T>(UPtrList<T>::set(i, nullptr));
}


template<class T>
inline void Foam::PtrList<T>::transfer(PtrList<T>& list)
{
    (this->ptrs_).free(); // free old pointers
    UPtrList<T>::transfer(list);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline void Foam::PtrList<T>::operator=(PtrList<T>&& list)
{
    this->transfer(list);
}


// ************************************************************************* //
