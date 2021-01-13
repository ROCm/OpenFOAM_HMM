/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
#include "refPtr.H"
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
    UPtrList<T>(list)
{
    // Took ownership of the pointers
    list = reinterpret_cast<T*>(0);
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
template<class... Args>
inline void Foam::PtrList<T>::emplace_append(Args&&... args)
{
    UPtrList<T>::append(new T(std::forward<Args>(args)...));
}


template<class T>
inline void Foam::PtrList<T>::append(T* ptr)
{
    UPtrList<T>::append(ptr);
}


template<class T>
inline void Foam::PtrList<T>::append(autoPtr<T>& ptr)
{
    UPtrList<T>::append(ptr.release());
}


template<class T>
inline void Foam::PtrList<T>::append(autoPtr<T>&& ptr)
{
    UPtrList<T>::append(ptr.release());
}


template<class T>
inline void Foam::PtrList<T>::append(std::unique_ptr<T>&& ptr)
{
    UPtrList<T>::append(ptr.release());
}


template<class T>
inline void Foam::PtrList<T>::append(const refPtr<T>& ptr)
{
    UPtrList<T>::append(ptr.ptr());
}


template<class T>
inline void Foam::PtrList<T>::append(const tmp<T>& ptr)
{
    UPtrList<T>::append(ptr.ptr());
}


template<class T>
template<class... Args>
inline Foam::autoPtr<T> Foam::PtrList<T>::emplace
(
    const label i,
    Args&&... args
)
{
    return set(i, new T(std::forward<Args>(args)...));
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrList<T>::set(const label i, T* ptr)
{
    // UPtrList::set returns a nullptr if trying to set with the same
    // pointer (again). This prevents the autoPtr from managing the
    // memory (avoids possible double deletion).

    return autoPtr<T>(UPtrList<T>::set(i, ptr));
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrList<T>::set
(
    const label i,
    autoPtr<T>& ptr
)
{
    return set(i, ptr.release());
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrList<T>::set
(
    const label i,
    autoPtr<T>&& ptr
)
{
    return set(i, ptr.release());
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrList<T>::set
(
    const label i,
    std::unique_ptr<T>&& ptr
)
{
    return set(i, ptr.release());
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrList<T>::set
(
    const label i,
    const refPtr<T>& ptr
)
{
    return set(i, ptr.ptr());
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrList<T>::set
(
    const label i,
    const tmp<T>& ptr
)
{
    return set(i, ptr.ptr());
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
