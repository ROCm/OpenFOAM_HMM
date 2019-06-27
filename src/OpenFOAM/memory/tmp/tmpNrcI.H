/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2016-2017 OpenFOAM Foundation
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
inline Foam::tmpNrc<T> Foam::tmpNrc<T>::New(Args&&... args)
{
    return tmpNrc<T>(new T(std::forward<Args>(args)...));
}


template<class T>
template<class U, class... Args>
inline Foam::tmpNrc<T> Foam::tmpNrc<T>::NewFrom(Args&&... args)
{
    return tmpNrc<T>(new U(std::forward<Args>(args)...));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline constexpr Foam::tmpNrc<T>::tmpNrc() noexcept
:
    ptr_(nullptr),
    type_(PTR)
{}


template<class T>
inline constexpr Foam::tmpNrc<T>::tmpNrc(std::nullptr_t) noexcept
:
    ptr_(nullptr),
    type_(PTR)
{}


template<class T>
inline Foam::tmpNrc<T>::tmpNrc(T* p) noexcept
:
    ptr_(p),
    type_(PTR)
{}


template<class T>
inline Foam::tmpNrc<T>::tmpNrc(const T& obj) noexcept
:
    ptr_(const_cast<T*>(&obj)),
    type_(CREF)
{}


template<class T>
inline Foam::tmpNrc<T>::tmpNrc(tmpNrc<T>&& t)
:
    ptr_(t.ptr_),
    type_(t.type_)
{
    t.ptr_ = nullptr;
    t.type_ = PTR;
}


template<class T>
inline Foam::tmpNrc<T>::tmpNrc(const tmpNrc<T>& t)
:
    ptr_(t.ptr_),
    type_(t.type_)
{
    if (isTmp())
    {
        if (ptr_)
        {
            t.type_ = CREF;
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
inline Foam::tmpNrc<T>::tmpNrc(const tmpNrc<T>& t, bool reuse)
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
                t.type_ = CREF;
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
inline Foam::tmpNrc<T>::~tmpNrc()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline bool Foam::tmpNrc<T>::isTmp() const noexcept
{
    return type_ == PTR;
}


template<class T>
inline bool Foam::tmpNrc<T>::empty() const noexcept
{
    return (!ptr_ && isTmp());
}


template<class T>
inline bool Foam::tmpNrc<T>::valid() const noexcept
{
    return (ptr_ || type_ == CREF);
}


template<class T>
inline bool Foam::tmpNrc<T>::movable() const noexcept
{
    return (type_ == PTR && ptr_);
}


template<class T>
inline Foam::word Foam::tmpNrc<T>::typeName() const
{
    return "tmpNrc<" + word(typeid(T).name()) + '>';
}


template<class T>
inline T* Foam::tmpNrc<T>::get() noexcept
{
    return ptr_; // non-const pointer
}


template<class T>
inline const T* Foam::tmpNrc<T>::get() const noexcept
{
    return ptr_; // const pointer
}


template<class T>
inline const T& Foam::tmpNrc<T>::cref() const
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
inline T& Foam::tmpNrc<T>::ref() const
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
inline T& Foam::tmpNrc<T>::constCast() const
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
inline T* Foam::tmpNrc<T>::ptr() const
{
    if (isTmp())
    {
        if (!ptr_)
        {
            FatalErrorInFunction
                << typeName() << " deallocated"
                << abort(FatalError);
        }

        T* ptr = ptr_;
        ptr_ = nullptr;

        return ptr;
    }

    return ptr_->clone().ptr();
}


template<class T>
inline void Foam::tmpNrc<T>::clear() const noexcept
{
    if (isTmp() && ptr_)
    {
        delete ptr_;
        ptr_ = nullptr;
    }
}


template<class T>
inline void Foam::tmpNrc<T>::reset() noexcept
{
    clear();
    ptr_ = nullptr;
    type_ = PTR;
}


template<class T>
inline void Foam::tmpNrc<T>::reset(T* p) noexcept
{
    clear();
    ptr_ = p;
    type_ = PTR;
}


template<class T>
inline void Foam::tmpNrc<T>::cref(const T& obj) noexcept
{
    clear();
    ptr_ = const_cast<T*>(&obj);
    type_ = CREF;
}


template<class T>
inline void Foam::tmpNrc<T>::swap(tmpNrc<T>& other) noexcept
{
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
inline const T& Foam::tmpNrc<T>::operator()() const
{
    return cref();
}


template<class T>
inline Foam::tmpNrc<T>::operator const T&() const
{
    return cref();
}


template<class T>
inline const T* Foam::tmpNrc<T>::operator->() const
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
inline T* Foam::tmpNrc<T>::operator->()
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
inline void Foam::tmpNrc<T>::operator=(T* p)
{
    clear();

    if (!p)
    {
        FatalErrorInFunction
            << "Attempted copy of a deallocated " << typeName()
            << abort(FatalError);
    }

    ptr_ = p;
    type_ = PTR;
}


template<class T>
inline void Foam::tmpNrc<T>::operator=(const tmpNrc<T>& t)
{
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
inline Foam::tmpNrc<T>::operator tmp<T>()
{
    if (isTmp())
    {
        return tmp<T>(ptr());
    }
    else
    {
        return tmp<T>(cref());
    }
}


// ************************************************************************* //
