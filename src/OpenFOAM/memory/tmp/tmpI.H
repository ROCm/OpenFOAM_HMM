/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

// * * * * * * * * * * * * * Private Member Operators  * * * * * * * * * * * //

template<class T>
inline void Foam::tmp<T>::incrCount()
{
    ptr_->operator++();

    if (ptr_->count() > 1)
    {
        FatalErrorInFunction
            << "Attempt to create more than 2 tmp's referring to"
               " the same object of type "
            << tmp<T>::typeName()
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


template<class T>
inline Foam::word Foam::tmp<T>::typeName()
{
    return "tmp<" + word(typeid(T).name()) + '>';
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
            << "Attempted construction of a "
            << this->typeName()
            << " from non-unique pointer"
            << abort(FatalError);
    }
}


template<class T>
inline constexpr Foam::tmp<T>::tmp(const T& obj) noexcept
:
    ptr_(const_cast<T*>(&obj)),
    type_(CREF)
{}


template<class T>
inline Foam::tmp<T>::tmp(tmp<T>&& rhs) noexcept
:
    ptr_(rhs.ptr_),
    type_(rhs.type_)
{
    rhs.ptr_ = nullptr;
    rhs.type_ = PTR;
}


template<class T>
inline Foam::tmp<T>::tmp(const tmp<T>&& rhs) noexcept
:
    ptr_(rhs.ptr_),
    type_(rhs.type_)
{
    rhs.ptr_ = nullptr;
    rhs.type_ = PTR;
}


template<class T>
inline Foam::tmp<T>::tmp(const tmp<T>& rhs)
:
    ptr_(rhs.ptr_),
    type_(rhs.type_)
{
    if (type_ == PTR)
    {
        if (ptr_)
        {
            this->incrCount();
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
inline Foam::tmp<T>::tmp(const tmp<T>& rhs, bool reuse)
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
                this->incrCount();
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
inline Foam::tmp<T>::~tmp()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline bool Foam::tmp<T>::movable() const noexcept
{
    return (type_ == PTR && ptr_ && ptr_->unique());
}


template<class T>
inline const T& Foam::tmp<T>::cref() const
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
inline T& Foam::tmp<T>::ref() const
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
inline T& Foam::tmp<T>::constCast() const
{
    return const_cast<T&>(cref());
}


template<class T>
inline T* Foam::tmp<T>::ptr() const
{
    if (!ptr_)
    {
        FatalErrorInFunction
            << this->typeName() << " deallocated"
            << abort(FatalError);
    }

    if (type_ == PTR)
    {
        if (!ptr_->unique())
        {
            FatalErrorInFunction
                << "Attempt to acquire pointer to object referred to"
                << " by multiple temporaries of type "
                << this->typeName()
                << abort(FatalError);
        }

        // Release pointer
        T* p = ptr_;
        ptr_ = nullptr;

        return p;
    }

    return ptr_->clone().ptr();
}


template<class T>
inline void Foam::tmp<T>::clear() const noexcept
{
    if (type_ == PTR && ptr_)
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
inline void Foam::tmp<T>::cref(const tmp<T>& other) noexcept
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
inline void Foam::tmp<T>::cref(const T& obj) noexcept
{
    clear();
    ptr_ = const_cast<T*>(&obj);
    type_ = CREF;
}


template<class T>
inline void Foam::tmp<T>::cref(const T* p) noexcept
{
    clear();
    ptr_ = const_cast<T*>(p);
    type_ = (ptr_ ? CREF : PTR);
}


template<class T>
inline void Foam::tmp<T>::ref(T& obj) noexcept
{
    clear();
    ptr_ = &obj;
    type_ = REF;
}


template<class T>
inline void Foam::tmp<T>::ref(T* p) noexcept
{
    clear();
    ptr_ = p;
    type_ = (ptr_ ? REF : PTR);
}


template<class T>
inline void Foam::tmp<T>::swap(tmp<T>& other) noexcept
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
inline const T* Foam::tmp<T>::operator->() const
{
    if (type_ == PTR && !ptr_)
    {
        FatalErrorInFunction
            << this->typeName() << " deallocated"
            << abort(FatalError);
    }

    return ptr_;
}


template<class T>
inline T* Foam::tmp<T>::operator->()
{
    if (is_const())
    {
        FatalErrorInFunction
            << "Attempt to cast const object to non-const: "
            << this->typeName()
            << abort(FatalError);
    }
    else if (type_ == PTR && !ptr_)
    {
        FatalErrorInFunction
            << this->typeName() << " deallocated"
            << abort(FatalError);
    }

    return ptr_;
}


template<class T>
inline void Foam::tmp<T>::operator=(const tmp<T>& other)
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


template<class T>
inline void Foam::tmp<T>::operator=(T* p)
{
    if (!p)
    {
        FatalErrorInFunction
            << "Attempted copy of a deallocated "
            << this->typeName()
            << abort(FatalError);
    }
    else if (!p->unique())
    {
        FatalErrorInFunction
            << "Attempted assignment of a "
            << this->typeName()
            << " to non-unique pointer"
            << abort(FatalError);
    }

    reset(p);
}


template<class T>
inline void Foam::tmp<T>::operator=(std::nullptr_t) noexcept
{
    reset(nullptr);
}


// ************************************************************************* //
