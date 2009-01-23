/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA


\*---------------------------------------------------------------------------*/

#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<int nBits>
inline unsigned int Foam::PackedList<nBits>::max_bits()
{
    return sizeof(PackedStorage)*8 - 1;
}


template<int nBits>
inline unsigned int Foam::PackedList<nBits>::max_value()
{
    return ((1u << nBits) - 1);
}


template<int nBits>
inline unsigned int Foam::PackedList<nBits>::packing()
{
    return sizeof(PackedStorage)*8 / nBits;
}


template<int nBits>
inline Foam::label Foam::PackedList<nBits>::packedLength(const label nElem)
{
    return (nElem + packing() - 1) / packing();
}


template<int nBits>
inline void Foam::PackedList<nBits>::checkIndex(const label i) const
{
    if (!size_)
    {
        FatalErrorIn("PackedList<nBits>::checkIndex(const label)")
            << "attempt to access element from zero-sized list"
            << abort(FatalError);
    }
    else if (i < 0 || i >= size_)
    {
        FatalErrorIn("PackedList<nBits>::checkIndex(const label)")
            << "index " << i << " out of range 0 ... " << size_-1
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<int nBits>
inline Foam::PackedList<nBits>::PackedList()
:
    List<PackedStorage>(),
    size_(0)
{}


template<int nBits>
inline Foam::PackedList<nBits>::PackedList(const label size)
:
    List<PackedStorage>(packedLength(size), 0u),
    size_(size)
{}


template<int nBits>
inline Foam::PackedList<nBits>::PackedList(const PackedList<nBits>& lst)
:
    List<PackedStorage>(lst),
    size_(lst.size())
{}


template<int nBits>
inline Foam::PackedList<nBits>::PackedList(const Xfer<PackedList<nBits> >& lst)
{
    transfer(lst());
}


template<int nBits>
inline Foam::autoPtr<Foam::PackedList<nBits> >
Foam::PackedList<nBits>::clone() const
{
    return autoPtr<PackedList<nBits> >(new PackedList<nBits>(*this));
}


// * * * * * * * * * * * * * * * * Iterators * * * * * * * * * * * * * * * * //

template<int nBits>
inline Foam::PackedList<nBits>::iterator::iterator
(
    const PackedList<nBits>& lst,
    const label i
)
:
    list_(const_cast<PackedList<nBits>&>(lst)),
    elem_(i / packing()),
    offset_(i % packing())
{}


template<int nBits>
inline void Foam::PackedList<nBits>::iterator::operator=
(
    const iterator& iter
)
{
    elem_ = iter.elem_;
    offset_ = iter.offset_;
}


template<int nBits>
inline bool Foam::PackedList<nBits>::iterator::operator==
(
    const iterator& iter
) const
{
    return elem_ == iter.elem_ && offset_ == iter.offset_;
}


template<int nBits>
inline bool Foam::PackedList<nBits>::iterator::operator!=
(
    const iterator& iter
) const
{
    return !(operator==(iter));
}


template<int nBits>
inline bool Foam::PackedList<nBits>::iterator::operator=
(
    const unsigned int val
)
{
    const unsigned int startBit = nBits * offset_;
    const unsigned int mask = max_value() << startBit;

    unsigned int& stored = chunk();
    const unsigned int old = stored;

    stored &= ~mask;
    stored |= (val << startBit) & mask;

    return (old != stored);
}


template<int nBits>
inline Foam::label Foam::PackedList<nBits>::iterator::position() const
{
    return (elem_ * packing() + offset_);
}


template<int nBits>
inline unsigned int& Foam::PackedList<nBits>::iterator::chunk() const
{
    return list_.List<PackedStorage>::operator[](elem_);
}


template<int nBits>
inline typename Foam::PackedList<nBits>::iterator&
Foam::PackedList<nBits>::iterator::operator++()
{
    if (position() < list_.size())
    {
        offset_++;
        if (offset_ >= packing())
        {
            elem_++;
            offset_ = 0;
        }
    }

    return *this;
}


template<int nBits>
inline typename Foam::PackedList<nBits>::iterator
Foam::PackedList<nBits>::iterator::operator++(int)
{
    iterator old = *this;
    ++*this;
    return old;
}


template<int nBits>
inline unsigned int
Foam::PackedList<nBits>::iterator::operator()() const
{
    return (chunk() >> (nBits * offset_)) & max_value();
}


template<int nBits>
inline Foam::PackedList<nBits>::iterator::operator
unsigned int () const
{
    return (chunk() >> (nBits * offset_)) & max_value();
}


template<int nBits>
inline Foam::PackedList<nBits>::iterator::operator
bool() const
{
    return !!(chunk() >> (nBits * offset_)) & max_value();
}


template<int nBits>
inline typename Foam::PackedList<nBits>::iterator
Foam::PackedList<nBits>::begin()
{
    return iterator(*this, 0);
}


template<int nBits>
inline typename Foam::PackedList<nBits>::iterator
Foam::PackedList<nBits>::end()
{
    return iterator(*this, size());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<int nBits>
inline Foam::label Foam::PackedList<nBits>::size() const
{
    return size_;
}


template<int nBits>
inline bool Foam::PackedList<nBits>::empty() const
{
    return !size_;
}


template<int nBits>
inline void Foam::PackedList<nBits>::resize
(
    const label nElem,
    const unsigned int& val
)
{
    reserve(nElem, val);
    size_ = nElem;
}


template<int nBits>
inline void Foam::PackedList<nBits>::setSize
(
    const label newSize,
    const unsigned int& val
)
{
    resize(newSize, val);
}



template<int nBits>
inline Foam::label Foam::PackedList<nBits>::capacity() const
{
    return packing() * List<PackedStorage>::size();
}


template<int nBits>
inline void Foam::PackedList<nBits>::setCapacity(const label nElem)
{
    List<PackedStorage>::setSize(packedLength(nElem), 0u);

    // truncates addressable section too
    if (size_ > nElem)
    {
#       ifdef DEBUGList
        // clear old values (purely cosmetics for string output)
        for (label i = nElem; i < size_; i++)
        {
            set(i, 0);
        }
#       endif
        size_ = nElem;
    }
    else
    {
        // fill new elements
        for (label i = size_; i < nElem; i++)
        {
            set(i, 0);
        }
    }
}


template<int nBits>
inline void Foam::PackedList<nBits>::reserve
(
    const label nElem,
    const unsigned int& val
)
{
    label len = packedLength(nElem);

    // need more capacity?
    if (len > List<PackedStorage>::size())
    {
        List<PackedStorage>::setSize(len, 0);
    }

    // fill new elements
    for (label i = size_; i < nElem; i++)
    {
        set(i, val);
    }
}


template<int nBits>
inline void Foam::PackedList<nBits>::clear()
{
    size_ = 0;
}


template<int nBits>
inline void Foam::PackedList<nBits>::clearStorage()
{
    List<PackedStorage>::clear();
    size_ = 0;
}


template<int nBits>
inline void Foam::PackedList<nBits>::shrink()
{
    label len = packedLength(size_);

    // we have unused space?
    if (len < List<PackedStorage>::size())
    {
        List<PackedStorage>::setSize(len);
    }
}

template<int nBits>
inline const Foam::List<unsigned int>&
Foam::PackedList<nBits>::storage() const
{
    return static_cast<const List<PackedStorage>&>(*this);
}


template<int nBits>
inline void Foam::PackedList<nBits>::transfer(PackedList<nBits>& lst)
{
    size_ = lst.size_;
    lst.size_ = 0;

    List<PackedStorage>::transfer(lst);
}


template<int nBits>
inline Foam::Xfer<Foam::PackedList<nBits> >
Foam::PackedList<nBits>::xfer()
{
    return xferMove(*this);
}


template<int nBits>
inline unsigned int Foam::PackedList<nBits>::get(const label i) const
{
#   ifdef DEBUGList
    checkIndex(i);
#   endif

    return iterator(*this, i)();
}


template<int nBits>
inline unsigned int Foam::PackedList<nBits>::operator[](const label i) const
{
    if (i >= size_)
    {
        setSize(i + 1);
        return 0;
    }
    else
    {
        return iterator(*this, i)();
    }
}


template<int nBits>
inline bool Foam::PackedList<nBits>::set
(
    const label i,
    const unsigned int val
)
{
#   ifdef DEBUGList
    checkIndex(i);

    if (val & ~max_value())
    {
        FatalErrorIn("PackedList<T>::set(const label, const unsigned int)")
            << "value " << label(val)
            << " out-of-range 0 ... " << label(max_value())
            << " representable by " << nBits << " bits"
            << abort(FatalError);
    }
#   endif

    return (iterator(*this, i) = val);
}


template<int nBits>
inline void Foam::PackedList<nBits>::append(const unsigned int val)
{
    label elemI = size_;
    reserve(elemI + 1);
    size_++;

    iterator(*this, elemI) = val;
}


template<int nBits>
inline typename Foam::PackedList<nBits>::iterator
Foam::PackedList<nBits>::operator[](const label i)
{
    if (i >= size_)
    {
        setSize(i + 1);
    }

    return iterator(*this, i);
}


template<int nBits>
inline void Foam::PackedList<nBits>::operator=(const unsigned int val)
{
    if (val)
    {
        forAll(*this, elemI)
        {
            set(elemI, val);
        }
    }
    else
    {
        List<PackedStorage>::operator=(val);
    }
}


// ************************************************************************* //
