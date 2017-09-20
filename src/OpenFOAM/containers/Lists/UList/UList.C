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


template<class T>
Foam::labelRange Foam::UList<T>::validateRange
(
    std::initializer_list<label> start_size_pair
) const
{
    if (start_size_pair.size() != 2)
    {
        FatalErrorInFunction
            << "range specified with " << start_size_pair.size()
            << " elements instead of 2"
            << abort(FatalError);
    }

    auto iter = start_size_pair.begin();

    const label beg = *(iter++);
    const label sz  = *iter;

    return this->validateRange(labelRange(beg, sz));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::UList<T>::deepCopy(const UList<T>& a)
{
    if (a.size_ != this->size_)
    {
        FatalErrorInFunction
            << "ULists have different sizes: "
            << this->size_ << " " << a.size_
            << abort(FatalError);
    }

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
                List_ELEM((*this), vp, i) = List_ELEM(a, ap, i);
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
Foam::UList<T> Foam::UList<T>::operator[]
(
    std::initializer_list<label> start_size_pair
)
{
    const labelRange slice = validateRange(start_size_pair);

    return UList<T>(&(this->v_[slice.start()]), slice.size()); // SubList
}


template<class T>
const Foam::UList<T> Foam::UList<T>::operator[]
(
    std::initializer_list<label> start_size_range
) const
{
    // Restricted range
    const labelRange slice = validateRange(start_size_range);

    return UList<T>(&(this->v_[slice.start()]), slice.size()); // SubList
}


template<class T>
void Foam::UList<T>::operator=(const T& t)
{
    List_ACCESS(T, (*this), vp);
    List_FOR_ALL((*this), i)
    {
        List_ELEM((*this), vp, i) = t;
    }
}


template<class T>
void Foam::UList<T>::operator=(const zero)
{
    List_ACCESS(T, (*this), vp);
    List_FOR_ALL((*this), i)
    {
        List_ELEM((*this), vp, i) = Zero;
    }
}


// * * * * * * * * * * * * * * STL Member Functions  * * * * * * * * * * * * //

template<class T>
void Foam::UList<T>::swap(UList<T>& a)
{
    Swap(size_, a.size_);
    Swap(v_, a.v_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
std::streamsize Foam::UList<T>::byteSize() const
{
    if (!contiguous<T>())
    {
        FatalErrorInFunction
            << "Cannot return the binary size of a list of "
               "non-primitive elements"
            << abort(FatalError);
    }

    return this->size_*sizeof(T);
}


template<class T>
void Foam::sort(UList<T>& a)
{
    std::sort(a.begin(), a.end());
}


template<class T, class Cmp>
void Foam::sort(UList<T>& a, const Cmp& cmp)
{
    std::sort(a.begin(), a.end(), cmp);
}


template<class T>
void Foam::stableSort(UList<T>& a)
{
    std::stable_sort(a.begin(), a.end());
}


template<class T, class Cmp>
void Foam::stableSort(UList<T>& a, const Cmp& cmp)
{
    std::stable_sort(a.begin(), a.end(), cmp);
}


template<class T>
void Foam::shuffle(UList<T>& a)
{
    std::random_shuffle(a.begin(), a.end());
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
bool Foam::UList<T>::operator==(const UList<T>& a) const
{
    bool equal = (this->size_ == a.size_);
    if (!equal)
    {
        return false;
    }

    List_CONST_ACCESS(T, (*this), vp);
    List_CONST_ACCESS(T, (a), ap);

    List_FOR_ALL((*this), i)
    {
        equal = (List_ELEM((*this), vp, i) == List_ELEM((a), ap, i));
        if (!equal) break;
    }

    return equal;
}


template<class T>
bool Foam::UList<T>::operator!=(const UList<T>& a) const
{
    return !operator==(a);
}


template<class T>
bool Foam::UList<T>::operator<(const UList<T>& a) const
{
    for
    (
        const_iterator vi = begin(), ai = a.begin();
        vi < end() && ai < a.end();
        ++vi, ++ai
    )
    {
        if (*vi < *ai)
        {
            return true;
        }
        else if (*vi > *ai)
        {
            return false;
        }
    }

    // Contents look to be identical, or lists have different sizes
    return (this->size_ < a.size_);
}


template<class T>
bool Foam::UList<T>::operator>(const UList<T>& a) const
{
    return a.operator<(*this);
}


template<class T>
bool Foam::UList<T>::operator<=(const UList<T>& a) const
{
    return !operator>(a);
}


template<class T>
bool Foam::UList<T>::operator>=(const UList<T>& a) const
{
    return !operator<(a);
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "UListIO.C"

// ************************************************************************* //
