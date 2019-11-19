/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "error.H"
#include <typeinfo>

// * * * * * * * * * * * * * Private Member Operators  * * * * * * * * * * * //

template<class T>
inline void Foam::tmp<T>::operator++()
{
    ptr_->operator++();

    if (ptr_->count() > 1)
    {
        FatalErrorInFunction
            << "Attempt to create more than 2 tmp's referring to"
               " the same object of type " << typeName()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class T>
template<class... Args>
inline Foam::tmp<T> Foam::tmp<T>::New(Args&&... args)
{
    return tmp<T>(new T(std::forward<Args>(args)...));
}


template<class T>
template<class U, class... Args>
inline Foam::tmp<T> Foam::tmp<T>::NewFrom(Args&&... args)
{
    return tmp<T>(new U(std::forward<Args>(args)...));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline constexpr Foam::tmp<T>::tmp() noexcept
:
    ptr_(nullptr),
    type_(PTR)
{}


template<class T>
inline constexpr Foam::tmp<T>::tmp(std::nullptr_t) noexcept
:
    ptr_(nullptr),
    type_(PTR)
{}


template<class T>
inline Foam::tmp<T>::tmp(T* p)
:
    ptr_(p),
    type_(PTR)
{
    if (p && !p->unique())
    {
        FatalErrorInFunction
            << "Attempted construction of a " << typeName()
            << " from non-unique pointer"
            << abort(FatalError);
    }
}


template<class T>
inline Foam::tmp<T>::tmp(const T& obj) noexcept
:
    ptr_(const_cast<T*>(&obj)),
    type_(CREF)
{}


template<class T>
inline Foam::tmp<T>::tmp(tmp<T>&& t) noexcept
:
    ptr_(t.ptr_),
    type_(t.type_)
{
    t.ptr_ = nullptr;
    t.type_ = PTR;
}


template<class T>
inline Foam::tmp<T>::tmp(const tmp<T>&& t) noexcept
:
    ptr_(t.ptr_),
    type_(t.type_)
{
    t.ptr_ = nullptr;
    t.type_ = PTR;
}


template<class T>
inline Foam::tmp<T>::tmp(const tmp<T>& t)
:
    ptr_(t.ptr_),
    type_(t.type_)
{
    if (isTmp())
    {
        if (ptr_)
        {
            operator++();
        }
        else
        {
            FatalErrorInFunction
                << "Attempted copy of a deallocated " << typeName()
                << abort(FatalError);
        }
    }
}


template<class T>
inline Foam::tmp<T>::tmp(const tmp<T>& t, bool reuse)
:
    ptr_(t.ptr_),
    type_(t.type_)
{
    if (isTmp())
    {
        if (ptr_)
        {
            if (reuse)
            {
                t.ptr_ = nullptr; // t.type_ already set as PTR
            }
            else
            {
                operator++();
            }
        }
        else
        {
            FatalErrorInFunction
                << "Attempted copy of a deallocated " << typeName()
                << abort(FatalError);
        }
    }
}


template<class T>
inline Foam::tmp<T>::~tmp()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline bool Foam::tmp<T>::isTmp() const noexcept
{
    return type_ == PTR;
}


template<class T>
inline bool Foam::tmp<T>::empty() const noexcept
{
    return (!ptr_ && isTmp());
}


template<class T>
inline bool Foam::tmp<T>::valid() const noexcept
{
    return (ptr_ || type_ == CREF);
}


template<class T>
inline bool Foam::tmp<T>::movable() const noexcept
{
    return (type_ == PTR && ptr_ && ptr_->unique());
}


template<class T>
inline Foam::word Foam::tmp<T>::typeName() const
{
    return "tmp<" + word(typeid(T).name()) + '>';
}


template<class T>
inline T* Foam::tmp<T>::get() noexcept
{
    return ptr_; // non-const pointer
}


template<class T>
inline const T* Foam::tmp<T>::get() const noexcept
{
    return ptr_; // const pointer
}


template<class T>
inline const T& Foam::tmp<T>::cref() const
{
    if (isTmp())
    {
        if (!ptr_)
        {
            FatalErrorInFunction
                << typeName() << " deallocated"
                << abort(FatalError);
        }
    }

    return *ptr_; // const reference
}


template<class T>
inline T& Foam::tmp<T>::ref() const
{
    if (isTmp())
    {
        if (!ptr_)
        {
            FatalErrorInFunction
                << typeName() << " deallocated"
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Attempted non-const reference to const object from a "
            << typeName()
            << abort(FatalError);
    }

    return *ptr_; // non-const reference
}


template<class T>
inline T& Foam::tmp<T>::constCast() const
{
    if (isTmp() && !ptr_)
    {
        FatalErrorInFunction
            << typeName() << " deallocated"
            << abort(FatalError);
    }

    return const_cast<T&>(*ptr_);
}


template<class T>
inline T* Foam::tmp<T>::ptr() const
{
    if (isTmp())
    {
        if (!ptr_)
        {
            FatalErrorInFunction
                << typeName() << " deallocated"
                << abort(FatalError);
        }
        else if (!ptr_->unique())
        {
            FatalErrorInFunction
                << "Attempt to acquire pointer to object referred to"
                << " by multiple temporaries of type " << typeName()
                << abort(FatalError);
        }

        T* ptr = ptr_;
        ptr_ = nullptr;

        return ptr;
    }

    return ptr_->clone().ptr();
}


template<class T>
inline void Foam::tmp<T>::clear() const noexcept
{
    if (isTmp() && ptr_)
    {
        if (ptr_->unique())
        {
            delete ptr_;
        }
        else
        {
            ptr_->operator--();
        }
        ptr_ = nullptr;
    }
}


template<class T>
inline void Foam::tmp<T>::reset() noexcept
{
    clear();
    ptr_ = nullptr;
    type_ = PTR;
}


template<class T>
inline void Foam::tmp<T>::reset(T* p) noexcept
{
    clear();
    ptr_ = p;
    type_ = PTR;
}


template<class T>
inline void Foam::tmp<T>::reset(tmp<T>&& other) noexcept
{
    if (&other == this)
    {
        return;  // Self-assignment is a no-op
    }

    clear();
    ptr_ = other.ptr_;
    type_ = other.type_;

    other.ptr_ = nullptr;
    other.type_ = PTR;
}


template<class T>
inline void Foam::tmp<T>::cref(const T& obj) noexcept
{
    clear();
    ptr_ = const_cast<T*>(&obj);
    type_ = CREF;
}


template<class T>
inline void Foam::tmp<T>::swap(tmp<T>& other) noexcept
{
    if (&other == this)
    {
        return;  // Self-swap is a no-op
    }

    // Copy/assign for pointer types
    T* p = ptr_;
    ptr_ = other.ptr_;
    other.ptr_ = p;

    refType t = type_;
    type_ = other.type_;
    other.type_ = t;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline const T& Foam::tmp<T>::operator()() const
{
    return cref();
}


template<class T>
inline Foam::tmp<T>::operator const T&() const
{
    return cref();
}


template<class T>
inline const T* Foam::tmp<T>::operator->() const
{
    if (!ptr_ && isTmp())
    {
        FatalErrorInFunction
            << typeName() << " deallocated"
            << abort(FatalError);
    }

    return ptr_;
}


template<class T>
inline T* Foam::tmp<T>::operator->()
{
    if (isTmp())
    {
        if (!ptr_)
        {
            FatalErrorInFunction
                << typeName() << " deallocated"
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Attempt to cast const object to non-const for a " << typeName()
            << abort(FatalError);
    }

    return ptr_;
}


template<class T>
inline Foam::tmp<T>::operator bool() const noexcept
{
    return (ptr_ || type_ == CREF);
}


template<class T>
inline void Foam::tmp<T>::operator=(T* p)
{
    clear();

    if (!p)
    {
        FatalErrorInFunction
            << "Attempted copy of a deallocated " << typeName()
            << abort(FatalError);
    }
    else if (!p->unique())
    {
        FatalErrorInFunction
            << "Attempted assignment of a " << typeName()
            << " to non-unique pointer"
            << abort(FatalError);
    }

    ptr_ = p;
    type_ = PTR;
}


template<class T>
inline void Foam::tmp<T>::operator=(const tmp<T>& t)
{
    if (&t == this)
    {
        return;  // Self-assignment is a no-op
    }

    clear();

    if (t.isTmp())
    {
        ptr_ = t.ptr_;
        type_ = PTR;
        t.ptr_ = nullptr;

        if (!ptr_)
        {
            FatalErrorInFunction
                << "Attempted assignment to a deallocated " << typeName()
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Attempted assignment to a const reference to an object"
            << " of type " << typeid(T).name()
            << abort(FatalError);
    }
}


template<class T>
inline void Foam::tmp<T>::operator=(tmp<T>&& other) noexcept
{
    if (&other == this)
    {
        return;  // Self-assignment is a no-op
    }

    clear();
    ptr_ = other.ptr_;
    type_ = other.type_;

    other.ptr_ = nullptr;
    other.type_ = PTR;
}


// ************************************************************************* //
