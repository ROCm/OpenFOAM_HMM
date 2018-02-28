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

#include "nullObject.H"
#include <utility>

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

template<class T>
inline const Foam::Xfer<T>& Foam::Xfer<T>::null()
{
    return NullObjectRef<Xfer<T>>();
}


template<class T>
template<class... Args>
inline Foam::Xfer<T> Foam::Xfer<T>::New(Args&&... args)
{
    return Xfer<T>(new T(std::forward<Args>(args)...));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::Xfer<T>::Xfer(T* p)
:
    ptr_(p ? p : new T)
{}


template<class T>
inline Foam::Xfer<T>::Xfer(const Xfer<T>& xf)
:
    ptr_(new T)
{
    T* p = xf.get();

    if (p && Foam::notNull(p))
    {
        ptr_->transfer(*p);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline bool Foam::Xfer<T>::valid() const
{
    return (ptr_ && Foam::notNull(ptr_));
}


template<class T>
inline T* Foam::Xfer<T>::get() const
{
    return ptr_;
}


template<class T>
inline T&& Foam::Xfer<T>::rvalue() const
{
    // Need check for valid() ?
    return std::move(*ptr_);
}


template<class T>
inline void Foam::Xfer<T>::swap(Xfer<T>& other) noexcept
{
    // Swap pointers
    T* p = ptr_;
    ptr_ = other.ptr;
    other.ptr_ = p;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class T>
inline Foam::Xfer<T>::~Xfer()
{
    if (ptr_ && Foam::notNull(ptr_))
    {
        delete ptr_;
    }
    ptr_ = nullptr;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline void Foam::Xfer<T>::operator=(const Xfer<T>& xf)
{
    // Silently ignore attempted copy to self
    if (this != &xf)
    {
        ptr_->transfer(*xf);
    }
}


template<class T>
inline T& Foam::Xfer<T>::operator*() const
{
    return *ptr_;
}


template<class T>
inline T* Foam::Xfer<T>::operator->() const
{
    return ptr_;
}


template<class T>
inline T& Foam::Xfer<T>::operator()() const
{
    return *ptr_;
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class T>
inline Foam::Xfer<T> Foam::xferCopy(const T& obj)
{
    T* ptr = new T(obj);
    return Foam::Xfer<T>(ptr);
}


template<class T, class From>
inline Foam::Xfer<T> Foam::xferCopyTo(const From& obj)
{
    T* ptr = new T(obj);
    return Foam::Xfer<T>(ptr);
}


template<class T>
inline Foam::Xfer<T> Foam::xferMove(T& obj)
{
    T* ptr = new T;
    ptr->transfer(obj);
    return Foam::Xfer<T>(ptr);
}


template<class T, class From>
inline Foam::Xfer<T> Foam::xferMoveTo(From& obj)
{
    T* ptr = new T;
    ptr->transfer(obj);
    return Foam::Xfer<T>(ptr);
}


// ************************************************************************* //
