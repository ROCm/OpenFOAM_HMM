/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "error.H"
#include <typeinfo>

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class T>
template<class... Args>
inline Foam::autoPtr<T> Foam::autoPtr<T>::New(Args&&... args)
{
    return autoPtr<T>(new T(std::forward<Args>(args)...));
}


template<class T>
template<class U, class... Args>
inline Foam::autoPtr<T> Foam::autoPtr<T>::NewFrom(Args&&... args)
{
    return autoPtr<T>(new U(std::forward<Args>(args)...));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline constexpr Foam::autoPtr<T>::autoPtr() noexcept
:
    ptr_(nullptr)
{}


template<class T>
inline constexpr Foam::autoPtr<T>::autoPtr(std::nullptr_t) noexcept
:
    ptr_(nullptr)
{}


template<class T>
inline Foam::autoPtr<T>::autoPtr(T* p) noexcept
:
    ptr_(p)
{}


template<class T>
inline Foam::autoPtr<T>::autoPtr(autoPtr<T>&& ap) noexcept
:
    ptr_(ap.release())
{}


template<class T>
template<class U>
inline Foam::autoPtr<T>::autoPtr(autoPtr<U>&& ap)
:
    ptr_(ap.release())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class T>
inline Foam::autoPtr<T>::~autoPtr() noexcept
{
    reset(nullptr);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline T* Foam::autoPtr<T>::release() noexcept
{
    T* p = ptr_;
    ptr_ = nullptr;
    return p;
}


template<class T>
inline void Foam::autoPtr<T>::reset(T* p) noexcept
{
    if (ptr_) delete ptr_;
    ptr_ = p;
}


template<class T>
inline void Foam::autoPtr<T>::reset(autoPtr<T>&& other) noexcept
{
    if (&other != this)
    {
        // Ignore self-assignment
        reset(other.release());
    }
}


template<class T>
inline void Foam::autoPtr<T>::swap(autoPtr<T>& other) noexcept
{
    // Swap is just copy/assign for pointer and enum types
    // Self-swap is effectively ignored
    T* p = ptr_;
    ptr_ = other.ptr_;
    other.ptr_ = p;
}


template<class T>
template<class... Args>
inline Foam::autoPtr<T> Foam::autoPtr<T>::clone(Args&&... args) const
{
    if (ptr_)
    {
        return autoPtr<T>(ptr_->clone(std::forward<Args>(args)...).ptr());
    }

    return autoPtr<T>();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline T& Foam::autoPtr<T>::operator*()
{
    if (!ptr_)
    {
        FatalErrorInFunction
            << "unallocated autoPtr of type " << typeid(T).name()
            << abort(FatalError);
    }
    return *ptr_;
}


template<class T>
inline const T& Foam::autoPtr<T>::operator*() const
{
    return const_cast<autoPtr<T>*>(this)->operator*();
}


template<class T>
inline T* Foam::autoPtr<T>::operator->()
{
    if (!ptr_)
    {
        FatalErrorInFunction
            << "unallocated autoPtr of type " << typeid(T).name()
            << abort(FatalError);
    }
    return ptr_;
}


template<class T>
inline const T* Foam::autoPtr<T>::operator->() const
{
    return const_cast<autoPtr<T>*>(this)->operator->();
}


template<class T>
inline T& Foam::autoPtr<T>::operator()()
{
    return operator*();
}


template<class T>
inline const T& Foam::autoPtr<T>::operator()() const
{
    return operator*();
}


template<class T>
inline void Foam::autoPtr<T>::operator=(autoPtr<T>&& ap) noexcept
{
    if (this != &ap)
    {
        // Ignore self-assignment
        reset(ap.release());
    }
}


template<class T>
template<class U>
inline void Foam::autoPtr<T>::operator=(autoPtr<U>&& ap) noexcept
{
    if (this != &ap)
    {
        // Ignore self-assignment
        reset(ap.release());
    }
}


template<class T>
inline void Foam::autoPtr<T>::operator=(std::nullptr_t) noexcept
{
    reset(nullptr);
}


// ************************************************************************* //
