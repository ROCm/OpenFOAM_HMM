/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
inline Foam::refPtr<T> Foam::refPtr<T>::New(Args&&... args)
{
    return refPtr<T>(new T(std::forward<Args>(args)...));
}


template<class T>
template<class U, class... Args>
inline Foam::refPtr<T> Foam::refPtr<T>::NewFrom(Args&&... args)
{
    return refPtr<T>(new U(std::forward<Args>(args)...));
}


template<class T>
inline Foam::word Foam::refPtr<T>::typeName()
{
    return "refPtr<" + word(typeid(T).name()) + '>';
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline constexpr Foam::refPtr<T>::refPtr() noexcept
:
    ptr_(nullptr),
    type_(PTR)
{}


template<class T>
inline constexpr Foam::refPtr<T>::refPtr(std::nullptr_t) noexcept
:
    ptr_(nullptr),
    type_(PTR)
{}


template<class T>
inline constexpr Foam::refPtr<T>::refPtr(T* p) noexcept
:
    ptr_(p),
    type_(PTR)
{}


template<class T>
inline Foam::refPtr<T>::refPtr(autoPtr<T>&& rhs) noexcept
:
    refPtr<T>(rhs.release())
{}


template<class T>
inline Foam::refPtr<T>::refPtr(std::unique_ptr<T>&& rhs) noexcept
:
    refPtr<T>(rhs.release())
{}


template<class T>
inline constexpr Foam::refPtr<T>::refPtr(const T& obj) noexcept
:
    ptr_(const_cast<T*>(&obj)),
    type_(CREF)
{}


template<class T>
inline Foam::refPtr<T>::refPtr(refPtr<T>&& rhs) noexcept
:
    ptr_(rhs.ptr_),
    type_(rhs.type_)
{
    rhs.ptr_ = nullptr;
    rhs.type_ = PTR;
}


template<class T>
inline Foam::refPtr<T>::refPtr(const refPtr<T>& rhs)
:
    ptr_(rhs.ptr_),
    type_(rhs.type_)
{
    if (type_ == PTR)
    {
        if (ptr_)
        {
            rhs.type_ = REF;  // (shallow copy)
        }
        else
        {
            FatalErrorInFunction
                << "Attempted copy of a deallocated "
                << this->typeName()
                << abort(FatalError);
        }
    }
}


template<class T>
inline Foam::refPtr<T>::refPtr(const refPtr<T>& rhs, bool reuse)
:
    ptr_(rhs.ptr_),
    type_(rhs.type_)
{
    if (type_ == PTR)
    {
        if (ptr_)
        {
            if (reuse)
            {
                rhs.ptr_ = nullptr;
                // Note: rhs.type_ already set as PTR
            }
            else
            {
                rhs.type_ = REF;  // (shallow copy)
            }
        }
        else
        {
            FatalErrorInFunction
                << "Attempted copy of a deallocated "
                << this->typeName()
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class T>
inline Foam::refPtr<T>::~refPtr()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline bool Foam::refPtr<T>::movable() const noexcept
{
    return (type_ == PTR && ptr_);
}


template<class T>
inline const T& Foam::refPtr<T>::cref() const
{
    if (type_ == PTR && !ptr_)
    {
        FatalErrorInFunction
            << this->typeName() << " deallocated"
            << abort(FatalError);
    }

    return *ptr_;  // const reference
}


template<class T>
inline T& Foam::refPtr<T>::ref() const
{
    if (is_const())
    {
        FatalErrorInFunction
            << "Attempted non-const reference to const object: "
            << this->typeName()
            << abort(FatalError);
    }
    else if (type_ == PTR && !ptr_)
    {
        FatalErrorInFunction
            << this->typeName() << " deallocated"
            << abort(FatalError);
    }

    return *ptr_;  // non-const reference
}


template<class T>
inline T& Foam::refPtr<T>::constCast() const
{
    return const_cast<T&>(cref());
}


template<class T>
inline Foam::refPtr<T> Foam::refPtr<T>::shallowClone() const noexcept
{
    refPtr<T> dup;

    if (ptr_)
    {
        dup.ptr_ = ptr_;
        dup.type_ = (is_const() ? CREF : REF);
    }

    return dup;
}


template<class T>
inline T* Foam::refPtr<T>::release() noexcept
{
    if (type_ == PTR)
    {
        T* p = ptr_;
        ptr_ = nullptr;
        return p;
    }

    return nullptr;
}


template<class T>
inline T* Foam::refPtr<T>::ptr() const
{
    if (!ptr_)
    {
        FatalErrorInFunction
            << this->typeName() << " deallocated"
            << abort(FatalError);
    }

    if (type_ == PTR)
    {
        // Release pointer
        T* p = ptr_;
        ptr_ = nullptr;

        return p;
    }

    return ptr_->clone().ptr();
}


template<class T>
inline void Foam::refPtr<T>::clear() const noexcept
{
    if (type_ == PTR && ptr_)
    {
        delete ptr_;
        ptr_ = nullptr;
    }
}


template<class T>
inline void Foam::refPtr<T>::reset(T* p) noexcept
{
    clear();
    ptr_ = p;
    type_ = PTR;
}


template<class T>
inline void Foam::refPtr<T>::reset(refPtr<T>&& other) noexcept
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
inline void Foam::refPtr<T>::cref(const refPtr<T>& other) noexcept
{
    if (&other == this)
    {
        return;  // Self-assignment is a no-op
    }

    clear();
    ptr_ = other.ptr_;
    type_ = (ptr_ ? CREF : PTR);
}


template<class T>
inline void Foam::refPtr<T>::cref(const T& obj) noexcept
{
    clear();
    ptr_ = const_cast<T*>(&obj);
    type_ = CREF;
}


template<class T>
inline void Foam::refPtr<T>::cref(const T* p) noexcept
{
    clear();
    ptr_ = const_cast<T*>(p);
    type_ = (ptr_ ? CREF : PTR);
}


template<class T>
inline void Foam::refPtr<T>::ref(T& obj) noexcept
{
    clear();
    ptr_ = &obj;
    type_ = REF;
}


template<class T>
inline void Foam::refPtr<T>::ref(T* p) noexcept
{
    clear();
    ptr_ = p;
    type_ = (ptr_ ? REF : PTR);
}


template<class T>
inline void Foam::refPtr<T>::swap(refPtr<T>& other) noexcept
{
    // Swap is just copy/assign for pointer and enum types
    // Self-swap is effectively ignored
    T* p = ptr_;
    ptr_ = other.ptr_;
    other.ptr_ = p;

    refType t = type_;
    type_ = other.type_;
    other.type_ = t;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline const T& Foam::refPtr<T>::operator*() const
{
    if (!ptr_)
    {
        FatalErrorInFunction
            << this->typeName() << " deallocated"
            << abort(FatalError);
    }
    return *ptr_;
}


template<class T>
inline T& Foam::refPtr<T>::operator*()
{
    if (is_const())
    {
        FatalErrorInFunction
            << "Attempt to cast const object to non-const: "
            << this->typeName()
            << abort(FatalError);
    }
    else if (!ptr_)
    {
        FatalErrorInFunction
            << this->typeName() << " deallocated"
            << abort(FatalError);
    }

    return *ptr_;
}


template<class T>
inline const T* Foam::refPtr<T>::operator->() const
{
    if (!ptr_)
    {
        FatalErrorInFunction
            << this->typeName() << " deallocated"
            << abort(FatalError);
    }

    return ptr_;
}


template<class T>
inline T* Foam::refPtr<T>::operator->()
{
    if (is_const())
    {
        FatalErrorInFunction
            << "Attempt to cast const object to non-const: "
            << this->typeName()
            << abort(FatalError);
    }
    else if (!ptr_)
    {
        FatalErrorInFunction
            << this->typeName() << " deallocated"
            << abort(FatalError);
    }

    return ptr_;
}


template<class T>
inline void Foam::refPtr<T>::operator=(const refPtr<T>& other)
{
    if (&other == this)
    {
        return;  // Self-assignment is a no-op
    }

    clear();

    if (other.type_ == PTR)
    {
        ptr_ = other.ptr_;
        type_ = PTR;
        other.ptr_ = nullptr;

        if (!ptr_)
        {
            FatalErrorInFunction
                << "Attempted assignment of a deallocated "
                << this->typeName()
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Attempted assignment of an object reference of type "
            << typeid(T).name()
            << abort(FatalError);
    }
}


template<class T>
inline void Foam::refPtr<T>::operator=(refPtr<T>&& other) noexcept
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
inline void Foam::refPtr<T>::operator=(T* p)
{
    if (!p)
    {
        FatalErrorInFunction
            << "Attempted copy of a deallocated "
            << this->typeName()
            << abort(FatalError);
    }

    reset(p);
}


template<class T>
inline void Foam::refPtr<T>::operator=(std::nullptr_t) noexcept
{
    reset(nullptr);
}


template<class T>
inline Foam::refPtr<T>::operator tmp<T>()
{
    if (type_ == PTR)
    {
        return tmp<T>(ptr());
    }
    else
    {
        return tmp<T>(cref());
    }
}


// ************************************************************************* //
