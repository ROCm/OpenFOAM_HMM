/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "List.H"
#include "ListLoopM.H"
#include "FixedList.H"
#include "PtrList.H"
#include "SLList.H"
#include "contiguous.H"
#include <utility>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T>
void Foam::List<T>::doResize(const label len)
{
    if (len == this->size_)
    {
        return;
    }

    if (len > 0)
    {
        // With sign-check to avoid spurious -Walloc-size-larger-than
        T* nv = new T[len];

        const label overlap = min(this->size_, len);

        if (overlap)
        {
            #ifdef USEMEMCPY
            if (is_contiguous<T>::value)
            {
                std::memcpy
                (
                    static_cast<void*>(nv), this->v_, overlap*sizeof(T)
                );
            }
            else
            #endif
            {
                List_ACCESS(T, *this, vp);
                for (label i = 0; i < overlap; ++i)
                {
                    nv[i] = std::move(vp[i]);
                }
            }
        }

        clear();
        this->size_ = len;
        this->v_ = nv;
    }
    else
    {
        // Or only #ifdef FULLDEBUG
        if (len < 0)
        {
            FatalErrorInFunction
                << "bad size " << len
                << abort(FatalError);
        }
        // #endif

        clear();
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
Foam::List<T>::List(const label len)
:
    UList<T>(nullptr, len)
{
    if (len < 0)
    {
        FatalErrorInFunction
            << "bad size " << len
            << abort(FatalError);
    }

    doAlloc();
}


template<class T>
Foam::List<T>::List(const label len, const T& val)
:
    UList<T>(nullptr, len)
{
    if (len < 0)
    {
        FatalErrorInFunction
            << "bad size " << len
            << abort(FatalError);
    }

    if (len)
    {
        doAlloc();

        List_ACCESS(T, (*this), vp);
        for (label i=0; i < len; ++i)
        {
            vp[i] = val;
        }
    }
}


template<class T>
Foam::List<T>::List(const label len, const Foam::zero)
:
    UList<T>(nullptr, len)
{
    if (len < 0)
    {
        FatalErrorInFunction
            << "bad size " << len
            << abort(FatalError);
    }

    if (len)
    {
        doAlloc();

        List_ACCESS(T, (*this), vp);
        for (label i=0; i < len; ++i)
        {
            vp[i] = Zero;
        }
    }
}


template<class T>
Foam::List<T>::List(const Foam::one, const T& val)
:
    UList<T>(new T[1], 1)
{
    this->v_[0] = val;
}


template<class T>
Foam::List<T>::List(const Foam::one, T&& val)
:
    UList<T>(new T[1], 1)
{
    this->v_[0] = std::move(val);
}


template<class T>
Foam::List<T>::List(const Foam::one, const Foam::zero)
:
    UList<T>(new T[1], 1)
{
    this->v_[0] = Zero;
}


template<class T>
Foam::List<T>::List(const UList<T>& a)
:
    UList<T>(nullptr, a.size_)
{
    const label len = this->size_;

    if (len)
    {
        doAlloc();

        #ifdef USEMEMCPY
        if (is_contiguous<T>::value)
        {
            std::memcpy
            (
                static_cast<void*>(this->v_), a.v_, this->size_bytes()
            );
        }
        else
        #endif
        {
            List_ACCESS(T, (*this), vp);
            List_CONST_ACCESS(T, a, ap);
            for (label i = 0; i < len; ++i)
            {
                vp[i] = ap[i];
            }
        }
    }
}


template<class T>
Foam::List<T>::List(const List<T>& a)
:
    UList<T>(nullptr, a.size_)
{
    const label len = this->size_;

    if (len)
    {
        doAlloc();

        #ifdef USEMEMCPY
        if (is_contiguous<T>::value)
        {
            std::memcpy
            (
                static_cast<void*>(this->v_), a.v_, this->size_bytes()
            );
        }
        else
        #endif
        {
            List_ACCESS(T, (*this), vp);
            List_CONST_ACCESS(T, a, ap);
            for (label i = 0; i < len; ++i)
            {
                vp[i] = ap[i];
            }
        }
    }
}


template<class T>
Foam::List<T>::List(List<T>& a, bool reuse)
:
    UList<T>(nullptr, a.size_)
{
    if (reuse)
    {
        // Steal content
        this->v_ = a.v_;
        a.v_ = nullptr;
        a.size_ = 0;
        return;
    }

    const label len = this->size_;

    if (len)
    {
        doAlloc();

        #ifdef USEMEMCPY
        if (is_contiguous<T>::value)
        {
            std::memcpy
            (
                static_cast<void*>(this->v_), a.v_, this->size_bytes()
            );
        }
        else
        #endif
        {
            List_ACCESS(T, (*this), vp);
            List_CONST_ACCESS(T, a, ap);
            for (label i = 0; i < len; ++i)
            {
                vp[i] = ap[i];
            }
        }
    }
}


template<class T>
Foam::List<T>::List(const UList<T>& list, const labelUList& indices)
:
    UList<T>(nullptr, indices.size())
{
    const label len = indices.size();

    if (len)
    {
        doAlloc();

        List_ACCESS(T, (*this), vp);

        for (label i=0; i < len; ++i)
        {
            vp[i] = list[indices[i]];
        }
    }
}


template<class T>
template<unsigned N>
Foam::List<T>::List
(
    const UList<T>& list,
    const FixedList<label,N>& indices
)
:
    UList<T>(nullptr, label(N))
{
    const label len = label(N);

    doAlloc();

    List_ACCESS(T, (*this), vp);

    for (label i=0; i < len; ++i)
    {
        vp[i] = list[indices[i]];
    }
}


template<class T>
template<unsigned N>
Foam::List<T>::List(const FixedList<T, N>& list)
:
    UList<T>(nullptr, label(N))
{
    doAlloc();
    copyList(list);
}


template<class T>
Foam::List<T>::List(const PtrList<T>& list)
:
    UList<T>(nullptr, list.size())
{
    doAlloc();
    copyList(list);
}


template<class T>
Foam::List<T>::List(const SLList<T>& list)
:
    List<T>(list.begin(), list.end(), list.size())
{}


template<class T>
template<class Addr>
Foam::List<T>::List(const IndirectListBase<T, Addr>& list)
:
    UList<T>(nullptr, list.size())
{
    doAlloc();
    copyList(list);
}


template<class T>
Foam::List<T>::List(std::initializer_list<T> list)
:
    List<T>(list.begin(), list.end(), list.size())
{}


template<class T>
Foam::List<T>::List(List<T>&& list)
:
    UList<T>()
{
    // Can use transfer or swap to manage content
    transfer(list);
}


template<class T>
template<int SizeMin>
Foam::List<T>::List(DynamicList<T, SizeMin>&& list)
:
    UList<T>()
{
    transfer(list);
}


template<class T>
Foam::List<T>::List(SortableList<T>&& list)
:
    UList<T>()
{
    transfer(list);
}


template<class T>
Foam::List<T>::List(SLList<T>&& list)
:
    UList<T>()
{
    operator=(std::move(list));
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class T>
Foam::List<T>::~List()
{
    if (this->v_)
    {
        delete[] this->v_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::List<T>::resize(const label len, const T& val)
{
    label idx = this->size_;
    this->doResize(len);

    List_ACCESS(T, *this, vp);
    while (idx < len)
    {
        vp[idx] = val;
        ++idx;
    }
}


template<class T>
void Foam::List<T>::transfer(List<T>& list)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    // Clear and swap - could also check for self assignment
    clear();
    this->size_ = list.size_;
    this->v_ = list.v_;

    list.size_ = 0;
    list.v_ = nullptr;
}


template<class T>
template<int SizeMin>
void Foam::List<T>::transfer(DynamicList<T, SizeMin>& list)
{
    // Shrink the allocated space to the number of elements used
    list.shrink();
    transfer(static_cast<List<T>&>(list));

    // Ensure DynamicList has proper capacity=0 too
    list.clearStorage();
}


template<class T>
void Foam::List<T>::transfer(SortableList<T>& list)
{
    // Shrink away the sort indices
    list.shrink();
    transfer(static_cast<List<T>&>(list));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
void Foam::List<T>::operator=(const UList<T>& a)
{
    if (this == &a)
    {
        return;  // Self-assignment is a no-op
    }

    reAlloc(a.size_);

    const label len = this->size_;

    if (len)
    {
        #ifdef USEMEMCPY
        if (is_contiguous<T>::value)
        {
            std::memcpy
            (
                static_cast<void*>(this->v_), a.v_, this->size_bytes()
            );
        }
        else
        #endif
        {
            List_ACCESS(T, (*this), vp);
            List_CONST_ACCESS(T, a, ap);
            for (label i = 0; i < len; ++i)
            {
                vp[i] = ap[i];
            }
        }
    }
}


template<class T>
void Foam::List<T>::operator=(const List<T>& list)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    operator=(static_cast<const UList<T>&>(list));
}


template<class T>
void Foam::List<T>::operator=(const SLList<T>& list)
{
    const label len = list.size();

    reAlloc(len);

    if (len)
    {
        T* iter = this->begin();

        for (const T& val : list)
        {
            *iter = val;
            ++iter;
        }
    }
}


template<class T>
template<unsigned N>
void Foam::List<T>::operator=(const FixedList<T, N>& list)
{
    reAlloc(static_cast<label>(N));

    T* iter = this->begin();

    for (const T& val : list)
    {
        *iter = val;
        ++iter;
    }
}


template<class T>
template<class Addr>
void Foam::List<T>::operator=(const IndirectListBase<T, Addr>& list)
{
    const label len = list.size();

    reAlloc(len);

    if (len)
    {
        List_ACCESS(T, (*this), vp);

        for (label i=0; i < len; ++i)
        {
            vp[i] = list[i];
        }
    }
}


template<class T>
void Foam::List<T>::operator=(std::initializer_list<T> list)
{
    const label len = list.size();

    reAlloc(len);

    if (len)
    {
        T* iter = this->begin();

        for (const T& val : list)
        {
            *iter = val;
            ++iter;
        }
    }
}


template<class T>
void Foam::List<T>::operator=(List<T>&& list)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    transfer(list);
}


template<class T>
template<int SizeMin>
void Foam::List<T>::operator=(DynamicList<T, SizeMin>&& list)
{
    transfer(list);
}


template<class T>
void Foam::List<T>::operator=(SortableList<T>&& list)
{
    transfer(list);
}


template<class T>
void Foam::List<T>::operator=(SLList<T>&& list)
{
    label len = list.size();

    reAlloc(len);

    T* iter = this->begin();

    while (len--)
    {
        *iter = std::move(list.removeHead());
        ++iter;
    }

    list.clear();
}


// ************************************************************************* //
