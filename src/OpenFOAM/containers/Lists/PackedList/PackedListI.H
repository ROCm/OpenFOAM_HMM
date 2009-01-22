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

#ifndef PackedListI_H
#define PackedListI_H

#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<int nBits>
inline unsigned int Foam::PackedBitRef<nBits>::max_bits()
{
    return sizeof(unsigned int)*8 - 1;
}


template<int nBits>
inline unsigned int Foam::PackedBitRef<nBits>::max_value()
{
    return ((1u << nBits) - 1);
}


template<int nBits>
inline unsigned int Foam::PackedBitRef<nBits>::packing()
{
    return sizeof(unsigned int)*8 / nBits;
}


// Calculate underlying list size
template<int nBits>
inline Foam::label Foam::PackedList<nBits>::storageSize(const label sz)
{
    return (sz + packing() - 1) / packing();
}


template<int nBits>
inline Foam::label Foam::PackedList<nBits>::location
(
    const label i,
    label& offset
)
{
    // the offset is the start bit within the storage element
    offset = nBits * (i % packing());
    return i / packing();
}


// Check value is representable in nBits
template<int nBits>
inline void Foam::PackedList<nBits>::checkValue(const unsigned int val)
{
    if (val > max_value())
    {
        FatalErrorIn("PackedList<T>::checkValue(const unsigned int)")
            << "value " << label(val) << " out of range 0 ... "
            << label(max_value())
            << " representable by " << nBits << " bits"
            << abort(FatalError);
    }
}


// Check index i is within valid range (0 ... size-1).
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

// Null constructor
template<int nBits>
inline Foam::PackedList<nBits>::PackedList()
:
    List<unsigned int>(),
    size_(0)
{}


// Construct with given size.
template<int nBits>
inline Foam::PackedList<nBits>::PackedList(const label size)
:
    List<unsigned int>(storageSize(size), 0u),
    size_(size)
{}


template<int nBits>
inline Foam::PackedList<nBits>::PackedList(const PackedList<nBits>& lst)
:
    List<unsigned int>(lst),
    size_(lst.size())
{}


template<int nBits>
inline Foam::PackedList<nBits>::PackedList(const Xfer<PackedList<nBits> >& lst)
{
    transfer(lst());
}



template<int nBits>
inline Foam::autoPtr<Foam::PackedList<nBits> > Foam::PackedList<nBits>::clone() const
{
    return autoPtr<PackedList<nBits> >(new PackedList<nBits>(*this));
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
    const label newSize,
    const unsigned int& val
)
{
    reserve(newSize, val);
    size_ = newSize;
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
    return packing() * List<unsigned int>::size();
}


template<int nBits>
inline void Foam::PackedList<nBits>::setCapacity(const label nElem)
{
    List<unsigned int>::setSize(storageSize(nElem), 0);

    for (label i = size_; i < nElem; i++)
    {
        set(i, 0);
    }

    size_ = nElem;
}


template<int nBits>
inline void Foam::PackedList<nBits>::reserve
(
    const label nElem,
    const unsigned int& val
)
{
    label rawSize = storageSize(nElem);

    // need more capacity?
    if (rawSize > List<unsigned int>::size())
    {
        List<unsigned int>::setSize(rawSize, 0);
    }

    // fill new elements, don't rely on List::setSize
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
    List<unsigned int>::clear();
    size_ = 0;
}


template<int nBits>
inline void Foam::PackedList<nBits>::shrink()
{
    label rawSize = storageSize(size_);

    // we have unused space?
    if (rawSize < List<unsigned int>::size())
    {
        List<unsigned int>::setSize(rawSize);
    }
}


// Get value at i
template<int nBits>
inline unsigned int Foam::PackedList<nBits>::get(const label i) const
{
#   ifdef DEBUGList
    checkIndex(i);
#   endif

    label startBit;
    const unsigned int& elem = List<unsigned int>::operator[]
    (
        location(i, startBit)
    );

    // PackedBitRef<nBits> iter(*this, elem, startBit);

    // return PackedBitRef<nBits>(*this, elem, startBit)();

    return (elem >> startBit) & max_value();
}


template<int nBits>
inline unsigned int Foam::PackedList<nBits>::operator[](const label i) const
{
    return get(i);
}


// Set value at i
template<int nBits>
inline bool Foam::PackedList<nBits>::set(const label i, const unsigned int val)
{
#   ifdef DEBUGList
    checkIndex(i);
    checkValue(val);
#   endif

    label startBit;
    unsigned int& elem = List<unsigned int>::operator[]
    (
        location(i, startBit)
    );

    unsigned int oldElem = elem;

    unsigned int shiftedMask = max_value() << startBit;
    unsigned int shiftedVal  = (val & max_value()) << startBit;

    elem = (elem & ~shiftedMask) | shiftedVal;

    return elem != oldElem;
}


template<int nBits>
inline bool Foam::PackedList<nBits>::append(const unsigned int val)
{
    label elemI = size_;
    reserve(elemI + 1);
    size_++;

    return set(elemI, val);
}


template<int nBits>
inline const Foam::List<unsigned int>& Foam::PackedList<nBits>::storage() const
{
    return static_cast<const List<unsigned int>&>(*this);
}


template<int nBits>
inline void Foam::PackedList<nBits>::transfer(PackedList<nBits>& lst)
{
    size_ = lst.size_;
    lst.size_ = 0;

    List<unsigned int>::transfer(lst);
}



template<int nBits>
inline Foam::Xfer<Foam::PackedList<nBits> >
Foam::PackedList<nBits>::xfer()
{
    return xferMove(*this);
}



template<int nBits>
inline typename Foam::PackedBitRef<nBits>
Foam::PackedList<nBits>::operator[](const label i)
{
    if (i >= size_)
    {
        setSize(i + 1);
    }

    label startBit;
    unsigned int& elem = List<unsigned int>::operator[]
    (
        location(i, startBit)
    );

    return ::Foam::PackedBitRef<nBits>
    (
        // *this,
        elem, 
        startBit
    );
}


// Set all to val
template<int nBits>
inline void Foam::PackedList<nBits>::operator=(const unsigned int val)
{
#   ifdef DEBUGList
    checkValue(val);
#   endif

    if (val)
    {
        forAll(*this, i)
        {
            set(i, val);
        }
    }
    else
    {
        List<unsigned int>::operator=(val);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
