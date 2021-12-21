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

#include "UList.H"
#include "ListLoopM.H"
#include "contiguous.H"
#include "labelRange.H"

#include <algorithm>
#include <random>

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class T>
Foam::labelRange
Foam::UList<T>::validateRange(const labelRange& requestedRange) const
{
    const labelRange range(requestedRange.subset0(this->size()));

    #ifdef FULLDEBUG
    this->checkStart(range.start());
    this->checkSize(range.start() + range.size());
    #endif

    return range;
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
            << "Lists have different sizes: "
            << len << " != " << list.size() << nl
            << abort(FatalError);
    }
    else if (len)
    {
        #ifdef USEMEMCPY
        if (is_contiguous<T>::value)
        {
            std::memcpy
            (
                static_cast<void*>(this->v_), list.v_, this->size_bytes()
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


template<class T>
template<class Addr>
void Foam::UList<T>::deepCopy(const IndirectListBase<T, Addr>& list)
{
    const label len = this->size_;

    if (len != list.size())
    {
        FatalErrorInFunction
            << "Lists have different sizes: "
            << len << " != " << list.size() << nl
            << abort(FatalError);
    }
    else if (len)
    {
        List_ACCESS(T, (*this), lhs);
        for (label i = 0; i < len; ++i)
        {
            lhs[i] = list[i];
        }
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
void Foam::UList<T>::operator=(const T& val)
{
    const label len = this->size();

    List_ACCESS(T, (*this), vp);

    for (label i=0; i < len; ++i)
    {
        vp[i] = val;
    }
}


template<class T>
void Foam::UList<T>::operator=(const Foam::zero)
{
    const label len = this->size();

    List_ACCESS(T, (*this), vp);

    for (label i=0; i < len; ++i)
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
            << "Invalid for non-contiguous data types"
            << abort(FatalError);
    }
    return this->size_bytes();
}


template<class T>
Foam::label Foam::UList<T>::find(const T& val, label pos) const
{
    const label len = this->size();

    if (pos >= 0 && len)
    {
        List_CONST_ACCESS(T, (*this), list);

        while (pos < len)
        {
            if (list[pos] == val)
            {
                return pos;
            }

            ++pos;
        }
    }

    return -1;
}


template<class T>
Foam::label Foam::UList<T>::rfind(const T& val, label pos) const
{
    // pos == -1 has same meaning as std::string::npos - search from end
    if (pos < 0 || pos >= this->size())
    {
        pos = this->size()-1;
    }

    List_CONST_ACCESS(T, (*this), list);

    while (pos >= 0)
    {
        if (list[pos] == val)
        {
            return pos;
        }

        --pos;
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
    std::shuffle(a.begin(), a.end(), std::default_random_engine());
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


// ************************************************************************* //
