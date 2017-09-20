/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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
#include "IndirectList.H"
#include "UIndirectList.H"
#include "BiIndirectList.H"
#include "contiguous.H"

#include <utility>

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
Foam::List<T>::List(const label s)
:
    UList<T>(nullptr, s)
{
    if (this->size_ < 0)
    {
        FatalErrorInFunction
            << "bad size " << this->size_
            << abort(FatalError);
    }

    alloc();
}


template<class T>
Foam::List<T>::List(const label s, const T& val)
:
    UList<T>(nullptr, s)
{
    if (this->size_ < 0)
    {
        FatalErrorInFunction
            << "bad size " << this->size_
            << abort(FatalError);
    }

    alloc();

    if (this->size_)
    {
        List_ACCESS(T, (*this), vp);
        List_FOR_ALL((*this), i)
        {
            vp[i] = val;
        }
    }
}


template<class T>
Foam::List<T>::List(const label s, const zero)
:
    UList<T>(nullptr, s)
{
    if (this->size_ < 0)
    {
        FatalErrorInFunction
            << "bad size " << this->size_
            << abort(FatalError);
    }

    alloc();

    if (this->size_)
    {
        List_ACCESS(T, (*this), vp);
        List_FOR_ALL((*this), i)
        {
            vp[i] = Zero;
        }
    }
}


template<class T>
Foam::List<T>::List(const List<T>& a)
:
    UList<T>(nullptr, a.size_)
{
    if (this->size_)
    {
        alloc();

        #ifdef USEMEMCPY
        if (contiguous<T>())
        {
            memcpy(this->v_, a.v_, this->byteSize());
        }
        else
        #endif
        {
            List_ACCESS(T, (*this), vp);
            List_CONST_ACCESS(T, a, ap);
            List_FOR_ALL((*this), i)
            {
                vp[i] = ap[i];
            }
        }
    }
}


template<class T>
template<class T2>
Foam::List<T>::List(const List<T2>& a)
:
    UList<T>(nullptr, a.size())
{
    if (this->size_)
    {
        alloc();

        List_ACCESS(T, (*this), vp);
        List_CONST_ACCESS(T2, a, ap);
        List_FOR_ALL((*this), i)
        {
            vp[i] = T(ap[i]);
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
        // swap content
        this->v_ = a.v_;
        a.v_ = nullptr;
        a.size_ = 0;
    }
    else if (this->size_)
    {
        alloc();

        #ifdef USEMEMCPY
        if (contiguous<T>())
        {
            memcpy(this->v_, a.v_, this->byteSize());
        }
        else
        #endif
        {
            List_ACCESS(T, (*this), vp);
            List_CONST_ACCESS(T, a, ap);
            List_FOR_ALL((*this), i)
            {
                vp[i] = ap[i];
            }
        }
    }
}


template<class T>
Foam::List<T>::List(const UList<T>& lst, const labelUList& mapAddressing)
:
    UList<T>(nullptr, mapAddressing.size())
{
    if (this->size_)
    {
        alloc();

        List_ACCESS(T, (*this), vp);

        const label len = (*this).size();
        for (label i=0; i < len; ++i)
        {
            vp[i] = lst[mapAddressing[i]];
        }
    }
}


template<class T>
template<class InputIterator>
Foam::List<T>::List(InputIterator begIter, InputIterator endIter)
:
    List<T>(begIter, endIter, std::distance(begIter, endIter))
{}


template<class T>
template<unsigned Size>
Foam::List<T>::List(const FixedList<T, Size>& lst)
:
    UList<T>(nullptr, Size)
{
    allocCopyList(lst);
}


template<class T>
Foam::List<T>::List(const PtrList<T>& lst)
:
    UList<T>(nullptr, lst.size())
{
    allocCopyList(lst);
}


template<class T>
Foam::List<T>::List(const SLList<T>& lst)
:
    List<T>(lst.begin(), lst.end(), lst.size())
{}


template<class T>
Foam::List<T>::List(const UIndirectList<T>& lst)
:
    UList<T>(nullptr, lst.size())
{
    allocCopyList(lst);
}


template<class T>
Foam::List<T>::List(const BiIndirectList<T>& lst)
:
    UList<T>(nullptr, lst.size())
{
    allocCopyList(lst);
}


template<class T>
Foam::List<T>::List(std::initializer_list<T> lst)
:
    List<T>(lst.begin(), lst.end(), lst.size())
{}


template<class T>
Foam::List<T>::List(const Xfer<List<T>>& lst)
{
    transfer(lst());
}


template<class T>
Foam::List<T>::List(List<T>&& lst)
:
    UList<T>(nullptr, 0)
{
    // Can use transfer or swap to manage content
    transfer(lst);
}


template<class T>
template<int SizeMin>
Foam::List<T>::List(DynamicList<T, SizeMin>&& lst)
:
    UList<T>(nullptr, 0)
{
    transfer(lst);
}


template<class T>
Foam::List<T>::List(SortableList<T>&& lst)
:
    UList<T>(nullptr, 0)
{
    transfer(lst);
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
void Foam::List<T>::setSize(const label newSize)
{
    if (newSize < 0)
    {
        FatalErrorInFunction
            << "bad size " << newSize
            << abort(FatalError);
    }

    if (newSize != this->size_)
    {
        if (newSize > 0)
        {
            T* nv = new T[label(newSize)];

            const label overlap = min(this->size_, newSize);

            if (overlap)
            {
                #ifdef USEMEMCPY
                if (contiguous<T>())
                {
                    memcpy(nv, this->v_, overlap*sizeof(T));
                }
                else
                #endif
                {
                    // No speedup observed for copy assignment on simple types
                    List_ACCESS(T, *this, vp);
                    for (label i = 0; i < overlap; ++i)
                    {
                        nv[i] = std::move(vp[i]);
                    }
                }
            }

            clear();
            this->size_ = newSize;
            this->v_ = nv;
        }
        else
        {
            clear();
        }
    }
}


template<class T>
void Foam::List<T>::setSize(const label newSize, const T& val)
{
    const label oldSize = label(this->size_);
    this->setSize(newSize);

    List_ACCESS(T, *this, vp);
    for (label i = oldSize; i < newSize; ++i)
    {
        vp[i] = val;
    }
}


template<class T>
void Foam::List<T>::transfer(List<T>& lst)
{
    // Clear and swap - could also check for self assignment
    clear();
    this->size_ = lst.size_;
    this->v_ = lst.v_;

    lst.size_ = 0;
    lst.v_ = nullptr;
}


template<class T>
template<int SizeMin>
void Foam::List<T>::transfer(DynamicList<T, SizeMin>& lst)
{
    // Shrink the allocated space to the number of elements used
    lst.shrink();
    transfer(static_cast<List<T>&>(lst));

    // Ensure DynamicList has proper capacity=0 too
    lst.clearStorage();
}


template<class T>
void Foam::List<T>::transfer(SortableList<T>& lst)
{
    // Shrink away the sort indices
    lst.shrink();
    transfer(static_cast<List<T>&>(lst));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
void Foam::List<T>::operator=(const UList<T>& a)
{
    reAlloc(a.size_);

    if (this->size_)
    {
        #ifdef USEMEMCPY
        if (contiguous<T>())
        {
            memcpy(this->v_, a.v_, this->byteSize());
        }
        else
        #endif
        {
            List_ACCESS(T, (*this), vp);
            List_CONST_ACCESS(T, a, ap);
            List_FOR_ALL((*this), i)
            {
                vp[i] = ap[i];
            }
        }
    }
}


template<class T>
void Foam::List<T>::operator=(const List<T>& lst)
{
    if (this == &lst)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    operator=(static_cast<const UList<T>&>(lst));
}


template<class T>
void Foam::List<T>::operator=(const SLList<T>& lst)
{
    reAlloc(lst.size());

    if (this->size_)
    {
        label i = 0;
        for (auto iter = lst.cbegin(); iter != lst.cend(); ++iter)
        {
            this->operator[](i++) = *iter;
        }
    }
}


template<class T>
void Foam::List<T>::operator=(const UIndirectList<T>& lst)
{
    reAlloc(lst.size());
    copyList(lst);
}


template<class T>
void Foam::List<T>::operator=(const BiIndirectList<T>& lst)
{
    reAlloc(lst.size());
    copyList(lst);
}


template<class T>
void Foam::List<T>::operator=(std::initializer_list<T> lst)
{
    reAlloc(lst.size());

    label i = 0;
    for (const auto& val : lst)
    {
        this->operator[](i++) = val;
    }
}


template<class T>
void Foam::List<T>::operator=(List<T>&& lst)
{
    if (this == &lst)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }

    transfer(lst);
}


template<class T>
template<int SizeMin>
void Foam::List<T>::operator=(DynamicList<T, SizeMin>&& lst)
{
    transfer(lst);
}


template<class T>
void Foam::List<T>::operator=(SortableList<T>&& lst)
{
    transfer(lst);
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "ListIO.C"

// ************************************************************************* //
