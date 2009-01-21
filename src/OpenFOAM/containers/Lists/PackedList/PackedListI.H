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
inline unsigned int Foam::PackedList<nBits>::max_bits()
{
    return sizeof(unsigned int)*8 - 1;
}


template<int nBits>
inline unsigned int Foam::PackedList<nBits>::max_value()
{
    return ((1u << nBits) - 1);
}


template<int nBits>
inline unsigned int Foam::PackedList<nBits>::packing()
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
    List<unsigned int>(0),
    size_(0)
{}


// Construct with given size.
template<int nBits>
inline Foam::PackedList<nBits>::PackedList(const label size)
:
    List<unsigned int>(storageSize(size), 0u),
    size_(size)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<int nBits>
inline void Foam::PackedList<nBits>::resize(const label newSize)
{
    setSize(newSize);
}


template<int nBits>
inline void Foam::PackedList<nBits>::resize
(
    const label newSize,
    const unsigned int& val
)
{
    setSize(newSize, val);
}


template<int nBits>
inline Foam::label Foam::PackedList<nBits>::capacity() const
{
    return packing() * List<unsigned int>::size();
}


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
inline Foam::List<unsigned int>& Foam::PackedList<nBits>::storage()
{
    return static_cast<List<unsigned int>&>(*this);
}


template<int nBits>
inline const Foam::List<unsigned int>& Foam::PackedList<nBits>::storage() const
{
    return static_cast<const List<unsigned int>&>(*this);
}


template<int nBits>
inline Foam::Xfer<Foam::PackedList<nBits> >
Foam::PackedList<nBits>::xfer()
{
    return xferMove(*this);
}



template<int nBits>
inline Foam::PackedBitRef
Foam::PackedList<nBits>::operator[](const label i)
{
#   ifdef DEBUGList
    checkIndex(i);
#   endif

    label startBit;
    unsigned int& elem = List<unsigned int>::operator[]
    (
        location(i, startBit)
    );

    return ::Foam::PackedBitRef(elem, startBit, max_value());
}


// Set all to val
template<int nBits>
inline void Foam::PackedList<nBits>::operator=(const unsigned int val)
{
#   ifdef DEBUGList
    checkValue(val);
#   endif

    if (val == 0)
    {
        List<unsigned int>::operator=(val);
    }
    else
    {
        for (label i = 0; i < size_; i++)
        {
            set(i, val);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
