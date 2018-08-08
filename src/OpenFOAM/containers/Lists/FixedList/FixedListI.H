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
#include "SLList.H"
#include <type_traits>
#include <utility>

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class T, unsigned Size>
inline bool Foam::FixedList<T, Size>::uniform() const
{
    if (Size > 1)
    {
        const T& val = first();

        for (unsigned i=1; i<Size; ++i)
        {
            if (val != (*this)[i])
            {
                return false;
            }
        }

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, unsigned Size>
inline Foam::FixedList<T, Size>::FixedList(const T& val)
{
    for (unsigned i=0; i<Size; ++i)
    {
        v_[i] = val;
    }
}


template<class T, unsigned Size>
inline Foam::FixedList<T, Size>::FixedList(const zero)
{
    for (unsigned i=0; i<Size; ++i)
    {
        v_[i] = Zero;
    }
}


template<class T, unsigned Size>
inline Foam::FixedList<T, Size>::FixedList(const T list[Size])
{
    for (unsigned i=0; i<Size; ++i)
    {
        v_[i] = list[i];
    }
}


template<class T, unsigned Size>
inline Foam::FixedList<T, Size>::FixedList(const FixedList<T, Size>& list)
{
    for (unsigned i=0; i<Size; ++i)
    {
        v_[i] = list.v_[i];
    }
}


template<class T, unsigned Size>
inline Foam::FixedList<T, Size>::FixedList(FixedList<T, Size>&& list)
{
    for (unsigned i=0; i<Size; ++i)
    {
        v_[i] = std::move(list.v_[i]);
    }
}


template<class T, unsigned Size>
template<class InputIterator>
inline Foam::FixedList<T, Size>::FixedList
(
    InputIterator begIter,
    InputIterator endIter
)
{
    checkSize(std::distance(begIter, endIter));

    InputIterator iter = begIter;
    for (unsigned i=0; i<Size; ++i)
    {
        v_[i] = *iter;
        ++iter;
    }
}


template<class T, unsigned Size>
inline Foam::FixedList<T, Size>::FixedList(std::initializer_list<T> list)
{
    checkSize(list.size());

    auto iter = list.begin();
    for (unsigned i=0; i<Size; ++i)
    {
        v_[i] = *iter;
        ++iter;
    }
}


template<class T, unsigned Size>
inline Foam::FixedList<T, Size>::FixedList(const UList<T>& list)
{
    checkSize(list.size());

    for (unsigned i=0; i<Size; ++i)
    {
        v_[i] = list[i];
    }
}


template<class T, unsigned Size>
inline Foam::FixedList<T, Size>::FixedList(const SLList<T>& list)
{
    checkSize(list.size());

    auto iter = list.begin();
    for (unsigned i=0; i<Size; ++i)
    {
        v_[i] = *iter;
        ++iter;
    }
}


template<class T, unsigned Size>
inline Foam::autoPtr<Foam::FixedList<T, Size>>
Foam::FixedList<T, Size>::clone() const
{
    return autoPtr<FixedList<T, Size>>::New(*this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, unsigned Size>
inline const Foam::FixedList<T, Size>& Foam::FixedList<T, Size>::null()
{
    return NullObjectRef<FixedList<T, Size>>();
}


template<class T, unsigned Size>
inline Foam::label Foam::FixedList<T, Size>::fcIndex(const label i) const
{
    return (i == Size-1 ? 0 : i+1);
}


template<class T, unsigned Size>
inline const T& Foam::FixedList<T, Size>::fcValue(const label i) const
{
    return this->operator[](this->fcIndex(i));
}


template<class T, unsigned Size>
inline T& Foam::FixedList<T, Size>::fcValue(const label i)
{
    return this->operator[](this->fcIndex(i));
}


template<class T, unsigned Size>
inline Foam::label Foam::FixedList<T, Size>::rcIndex(const label i) const
{
    return (i ? i-1 : Size-1);
}


template<class T, unsigned Size>
inline const T& Foam::FixedList<T, Size>::rcValue(const label i) const
{
    return this->operator[](this->rcIndex(i));
}


template<class T, unsigned Size>
inline T& Foam::FixedList<T, Size>::rcValue(const label i)
{
    return this->operator[](this->rcIndex(i));
}


template<class T, unsigned Size>
inline void Foam::FixedList<T, Size>::checkStart(const label start) const
{
    if (start < 0 || (start && unsigned(start) >= Size))
    {
        // Note: always accept start=0, even for zero-sized lists
        FatalErrorInFunction
            << "start " << start << " out of range [0," << Size << ")"
            << abort(FatalError);
    }
}


template<class T, unsigned Size>
inline void Foam::FixedList<T, Size>::checkSize(const label size) const
{
    if (unsigned(size) != Size)
    {
        FatalErrorInFunction
            << "size " << size << " != " << Size
            << abort(FatalError);
    }
}


template<class T, unsigned Size>
inline void Foam::FixedList<T, Size>::checkIndex(const label i) const
{
    if (i < 0 || unsigned(i) >= Size)
    {
        FatalErrorInFunction
            << "index " << i << " out of range [0," << Size << ")"
            << abort(FatalError);
    }
}


template<class T, unsigned Size>
inline bool Foam::FixedList<T, Size>::found
(
    const T& val,
    const label start
) const
{
    return (this->find(val, start) >= 0);
}


template<class T, unsigned Size>
inline void Foam::FixedList<T, Size>::resize(const label n)
{
    #ifdef FULLDEBUG
    checkSize(n);
    #endif
}

template<class T, unsigned Size>
inline void Foam::FixedList<T, Size>::setSize(const label n)
{
    #ifdef FULLDEBUG
    checkSize(n);
    #endif
}


template<class T, unsigned Size>
inline void Foam::FixedList<T, Size>::swap(FixedList<T, Size>& list)
{
    for (unsigned i=0; i<Size; ++i)
    {
        Foam::Swap(v_[i], list.v_[i]);
    }
}


template<class T, unsigned Size>
inline void Foam::FixedList<T, Size>::transfer(FixedList<T, Size>& list)
{
    for (unsigned i=0; i<Size; ++i)
    {
        v_[i] = std::move(list[i]);
    }
}


template<class T, unsigned Size>
inline const T*
Foam::FixedList<T, Size>::cdata() const
{
    return v_;
}


template<class T, unsigned Size>
inline T*
Foam::FixedList<T, Size>::data()
{
    return v_;
}


template<class T, unsigned Size>
inline T& Foam::FixedList<T, Size>::first()
{
    return v_[0];
}


template<class T, unsigned Size>
inline const T& Foam::FixedList<T, Size>::first() const
{
    return v_[0];
}


template<class T, unsigned Size>
inline T& Foam::FixedList<T, Size>::last()
{
    return v_[Size-1];
}


template<class T, unsigned Size>
inline const T& Foam::FixedList<T, Size>::last() const
{
    return v_[Size-1];
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, unsigned Size>
inline T& Foam::FixedList<T, Size>::operator[](const label i)
{
    #ifdef FULLDEBUG
    checkIndex(i);
    #endif
    return v_[i];
}


template<class T, unsigned Size>
inline const T& Foam::FixedList<T, Size>::operator[](const label i) const
{
    #ifdef FULLDEBUG
    checkIndex(i);
    #endif
    return v_[i];
}


template<class T, unsigned Size>
inline void Foam::FixedList<T, Size>::operator=(const T list[Size])
{
    for (unsigned i=0; i<Size; ++i)
    {
        v_[i] = list[i];
    }
}

template<class T, unsigned Size>
inline void Foam::FixedList<T, Size>::operator=(const UList<T>& list)
{
    checkSize(list.size());

    for (unsigned i=0; i<Size; ++i)
    {
        v_[i] = list[i];
    }
}

template<class T, unsigned Size>
inline void Foam::FixedList<T, Size>::operator=(const SLList<T>& list)
{
    checkSize(list.size());

    auto iter = list.begin();
    for (unsigned i=0; i<Size; ++i)
    {
        v_[i] = *iter;
        ++iter;
    }
}

template<class T, unsigned Size>
inline void Foam::FixedList<T, Size>::operator=(std::initializer_list<T> list)
{
    checkSize(list.size());

    auto iter = list.begin();
    for (unsigned i=0; i<Size; ++i)
    {
        v_[i] = *iter;
        ++iter;
    }
}

template<class T, unsigned Size>
inline void Foam::FixedList<T, Size>::operator=(const T& val)
{
    for (unsigned i=0; i<Size; ++i)
    {
        v_[i] = val;
    }
}

template<class T, unsigned Size>
inline void Foam::FixedList<T, Size>::operator=(const FixedList<T, Size>& list)
{
    for (unsigned i=0; i<Size; ++i)
    {
        v_[i] = list.v_[i];
    }
}

template<class T, unsigned Size>
inline void Foam::FixedList<T, Size>::operator=(FixedList<T, Size>&& list)
{
    // No significant speedup observed for copy assignment on simple types,
    // use move assignment for generality with more complex types
    for (unsigned i=0; i<Size; ++i)
    {
        v_[i] = std::move(list.v_[i]);
    }
}


// * * * * * * * * * * * * * * STL Member Functions  * * * * * * * * * * * * //

template<class T, unsigned Size>
inline typename Foam::FixedList<T, Size>::iterator
Foam::FixedList<T, Size>::begin()
{
    return v_;
}


template<class T, unsigned Size>
inline typename Foam::FixedList<T, Size>::const_iterator
Foam::FixedList<T, Size>::begin() const
{
    return v_;
}


template<class T, unsigned Size>
inline typename Foam::FixedList<T, Size>::const_iterator
Foam::FixedList<T, Size>::cbegin() const
{
    return v_;
}


template<class T, unsigned Size>
inline typename Foam::FixedList<T, Size>::iterator
Foam::FixedList<T, Size>::end()
{
    return (v_ + Size);
}


template<class T, unsigned Size>
inline typename Foam::FixedList<T, Size>::const_iterator
Foam::FixedList<T, Size>::end() const
{
    return (v_ + Size);
}


template<class T, unsigned Size>
inline typename Foam::FixedList<T, Size>::const_iterator
Foam::FixedList<T, Size>::cend() const
{
    return (v_ + Size);
}


template<class T, unsigned Size>
inline typename Foam::FixedList<T, Size>::reverse_iterator
Foam::FixedList<T, Size>::rbegin()
{
    return reverse_iterator(end());
}


template<class T, unsigned Size>
inline typename Foam::FixedList<T, Size>::const_reverse_iterator
Foam::FixedList<T, Size>::rbegin() const
{
    return const_reverse_iterator(end());
}


template<class T, unsigned Size>
inline typename Foam::FixedList<T, Size>::const_reverse_iterator
Foam::FixedList<T, Size>::crbegin() const
{
    return const_reverse_iterator(end());
}


template<class T, unsigned Size>
inline typename Foam::FixedList<T, Size>::reverse_iterator
Foam::FixedList<T, Size>::rend()
{
    return reverse_iterator(begin());
}


template<class T, unsigned Size>
inline typename Foam::FixedList<T, Size>::const_reverse_iterator
Foam::FixedList<T, Size>::rend() const
{
    return const_reverse_iterator(begin());
}


template<class T, unsigned Size>
inline typename Foam::FixedList<T, Size>::const_reverse_iterator
Foam::FixedList<T, Size>::crend() const
{
    return const_reverse_iterator(begin());
}


template<class T, unsigned Size>
inline Foam::label Foam::FixedList<T, Size>::size() const
{
    return Size;
}


template<class T, unsigned Size>
inline Foam::label Foam::FixedList<T, Size>::max_size() const
{
    return Size;
}


template<class T, unsigned Size>
inline bool Foam::FixedList<T, Size>::empty() const
{
    return false;
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class T, unsigned Size>
inline void Foam::Swap(FixedList<T, Size>& lhs, FixedList<T, Size>& rhs)
{
    lhs.swap(rhs);
}


// ************************************************************************* //
