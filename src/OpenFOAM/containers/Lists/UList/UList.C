/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
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

#include "UList.H"
#include "ListLoopM.H"
#include "contiguous.H"
#include "labelRange.H"

#include <algorithm>

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class T>
Foam::labelRange Foam::UList<T>::validateRange(const labelRange& range) const
{
    const labelRange slice = range.subset0(this->size());

    #ifdef FULLDEBUG
    this->checkStart(slice.start());
    this->checkSize(slice.start() + slice.size());
    #endif

    return slice;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::UList<T>::moveFirst(const label i)
{
    checkIndex(i);

    for (label lower = 0; lower < i; ++lower)
    {
        Foam::Swap(this->operator[](lower), this->operator[](i));
    }
}


template<class T>
void Foam::UList<T>::moveLast(const label i)
{
    checkIndex(i);

    for (label upper = size()-1; upper > i; --upper)
    {
        Foam::Swap(this->operator[](i), this->operator[](upper));
    }
}


template<class T>
void Foam::UList<T>::swapFirst(const label i)
{
    checkIndex(i);

    if (i > 0)
    {
        Foam::Swap(this->operator[](0), this->operator[](i));
    }
}


template<class T>
void Foam::UList<T>::swapLast(const label i)
{
    checkIndex(i);

    const label upper = size()-1;

    if (i < upper)
    {
        Foam::Swap(this->operator[](i), this->operator[](upper));
    }
}


template<class T>
void Foam::UList<T>::deepCopy(const UList<T>& list)
{
    const label len = this->size_;

    if (len != list.size_)
    {
        FatalErrorInFunction
            << "ULists have different sizes: "
            << len << " " << list.size_
            << abort(FatalError);
    }
    else if (len)
    {
        #ifdef USEMEMCPY
        if (is_contiguous<T>::value)
        {
            std::memcpy
            (
                static_cast<void*>(this->v_), list.v_, this->byteSize()
            );
        }
        else
        #endif
        {
            List_ACCESS(T, (*this), lhs);
            List_CONST_ACCESS(T, list, rhs);
            for (label i = 0; i < len; ++i)
            {
                lhs[i] = rhs[i];
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
Foam::UList<T> Foam::UList<T>::operator[](const labelRange& range)
{
    const labelRange slice = validateRange(range);

    return UList<T>(&(this->v_[slice.start()]), slice.size()); // SubList
}


template<class T>
const Foam::UList<T> Foam::UList<T>::operator[](const labelRange& range) const
{
    const labelRange slice = validateRange(range);

    return UList<T>(&(this->v_[slice.start()]), slice.size()); // SubList
}


template<class T>
void Foam::UList<T>::operator=(const T& val)
{
    List_ACCESS(T, (*this), vp);
    List_FOR_ALL((*this), i)
    {
        vp[i] = val;
    }
}


template<class T>
void Foam::UList<T>::operator=(const zero)
{
    List_ACCESS(T, (*this), vp);
    List_FOR_ALL((*this), i)
    {
        vp[i] = Zero;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
std::streamsize Foam::UList<T>::byteSize() const
{
    if (!is_contiguous<T>::value)
    {
        FatalErrorInFunction
            << "Cannot return binary size of a list with non-primitive elements"
            << abort(FatalError);
    }

    return this->size_*sizeof(T);
}


template<class T>
Foam::label Foam::UList<T>::find(const T& val, const label start) const
{
    const label len = this->size();

    if (start >= 0 && len)
    {
        List_CONST_ACCESS(T, (*this), vp);

        for (label i = start; i < len; ++i)
        {
            if (vp[i] == val)
            {
                return i;
            }
        }
    }

    return -1;
}


template<class T>
Foam::label Foam::UList<T>::rfind(const T& val, const label pos) const
{
    List_CONST_ACCESS(T, (*this), vp);

    const label len1 = (this->size()-1);

    // pos == -1 has same meaning as std::string::npos - search from end
    for (label i = ((pos >= 0 && pos < len1) ? pos : len1); i >= 0; --i)
    {
        if (vp[i] == val)
        {
            return i;
        }
    }

    return -1;
}


template<class T>
void Foam::sort(UList<T>& a)
{
    std::sort(a.begin(), a.end());
}


template<class T, class Compare>
void Foam::sort(UList<T>& a, const Compare& comp)
{
    std::sort(a.begin(), a.end(), comp);
}


template<class T>
void Foam::stableSort(UList<T>& a)
{
    std::stable_sort(a.begin(), a.end());
}


template<class T, class Compare>
void Foam::stableSort(UList<T>& a, const Compare& comp)
{
    std::stable_sort(a.begin(), a.end(), comp);
}


template<class T>
void Foam::shuffle(UList<T>& a)
{
    std::random_shuffle(a.begin(), a.end());
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
bool Foam::UList<T>::operator==(const UList<T>& list) const
{
    const label len = this->size_;
    if (len != list.size_)
    {
        return false;
    }

    bool equal = true;

    List_CONST_ACCESS(T, (*this), lhs);
    List_CONST_ACCESS(T, (list), rhs);

    for (label i = 0; i < len; ++i)
    {
        equal = (lhs[i] == rhs[i]);
        if (!equal) break;
    }

    return equal;
}


template<class T>
bool Foam::UList<T>::operator!=(const UList<T>& list) const
{
    return !operator==(list);
}


template<class T>
bool Foam::UList<T>::operator<(const UList<T>& list) const
{
    for
    (
        const_iterator lhs = begin(), rhs = list.begin();
        lhs < end() && rhs < list.end();
        ++lhs, ++rhs
    )
    {
        if (*lhs < *rhs)
        {
            return true;
        }
        else if (*rhs < *lhs)
        {
            return false;
        }
    }

    // Contents look to be identical, or lists have different sizes
    return (this->size_ < list.size_);
}


template<class T>
bool Foam::UList<T>::operator>(const UList<T>& list) const
{
    return list.operator<(*this);
}


template<class T>
bool Foam::UList<T>::operator<=(const UList<T>& list) const
{
    return !list.operator<(*this);
}


template<class T>
bool Foam::UList<T>::operator>=(const UList<T>& list) const
{
    return !operator<(list);
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "UListIO.C"

// ************************************************************************* //
