/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017-2018 OpenCFD Ltd.
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

#include <climits>
#include "error.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<unsigned nBits>
inline constexpr unsigned int Foam::PackedList<nBits>::max_bits()
{
    return sizeof(StorageType)*CHAR_BIT - 1;
}


template<unsigned nBits>
inline constexpr unsigned int Foam::PackedList<nBits>::max_value()
{
    return (1u << nBits) - 1;
}


template<unsigned nBits>
inline constexpr unsigned int Foam::PackedList<nBits>::packing()
{
    return sizeof(StorageType)*CHAR_BIT / nBits;
}


template<unsigned nBits>
inline constexpr unsigned int Foam::PackedList<nBits>::maskLower
(
    unsigned offset
)
{
    // Return (1u << (nBits * offset)) - 1;
    // The next one works more reliably with overflows
    // eg, when compiled without optimization
    return (~0u >> (sizeof(StorageType)*CHAR_BIT - nBits * offset));
}


template<unsigned nBits>
inline constexpr Foam::label Foam::PackedList<nBits>::packedLength
(
    const label nElem
)
{
    return (nElem + packing() - 1) / packing();
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{
    // Template specialization for bool entries
    template<> inline unsigned int PackedList<1>::readValue(Istream& is)
    {
        return readBool(is);
    }

    // Template specialization for bool entries
    template<> inline void PackedList<1>::setPair(Istream& is)
    {
        set(readLabel(is), true);
    }
}


template<unsigned nBits>
inline unsigned int Foam::PackedList<nBits>::readValue(Istream& is)
{
    const unsigned int val = readLabel(is);

    if (val > max_value())
    {
        FatalIOErrorInFunction(is)
            << "Out-of-range value " << val << " for PackedList<" << nBits
            << ">. Maximum permitted value is " << max_value() << "."
            << exit(FatalIOError);
    }

    return val;
}


template<unsigned nBits>
inline void Foam::PackedList<nBits>::setPair(Istream& is)
{
    is.readBegin("Tuple2<label,unsigned int>");

    const label ind = readLabel(is);
    const unsigned int val = readLabel(is);

    is.readEnd("Tuple2<label,unsigned int>");

    if (val > max_value())
    {
        FatalIOErrorInFunction(is)
            << "Out-of-range value " << val << " for PackedList<" << nBits
            << "> at index " << ind
            << ". Maximum permitted value is " << max_value() << "."
            << exit(FatalIOError);
    }

    set(ind, val);

    is.check(FUNCTION_NAME);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<unsigned nBits>
inline Foam::PackedList<nBits>::PackedList()
:
    PackedListCore(),
    StorageList(),
    size_(0)
{}


template<unsigned nBits>
inline Foam::PackedList<nBits>::PackedList(const label size)
:
    PackedListCore(),
    StorageList(packedLength(size), 0u),
    size_(size)
{}


template<unsigned nBits>
inline Foam::PackedList<nBits>::PackedList
(
    const label size,
    const unsigned int val
)
:
    PackedListCore(),
    StorageList(packedLength(size), 0u),
    size_(size)
{
    if (val)
    {
        operator=(val);
    }
}


template<unsigned nBits>
inline Foam::PackedList<nBits>::PackedList(Istream& is)
:
    PackedListCore(),
    StorageList(),
    size_(0)
{
    read(is);
}


template<unsigned nBits>
inline Foam::PackedList<nBits>::PackedList(const PackedList<nBits>& lst)
:
    PackedListCore(),
    StorageList(lst),
    size_(lst.size_)
{}


template<unsigned nBits>
inline Foam::PackedList<nBits>::PackedList(PackedList<nBits>&& lst)
:
    PackedListCore(),
    StorageList(),
    size_(0)
{
    transfer(lst);
}


template<unsigned nBits>
inline Foam::PackedList<nBits>::PackedList(const labelUList& lst)
:
    PackedListCore(),
    StorageList(packedLength(lst.size()), 0u),
    size_(lst.size())
{
    const label len = lst.size();
    for (label i = 0; i < len; ++i)
    {
        set(i, lst[i]);
    }
}


template<unsigned nBits>
inline Foam::PackedList<nBits>::PackedList(const labelUIndList& lst)
:
    PackedListCore(),
    StorageList(packedLength(lst.size()), 0u),
    size_(lst.size())
{
    const label len = lst.size();
    for (label i = 0; i < len; ++i)
    {
        set(i, lst[i]);
    }
}


template<unsigned nBits>
inline Foam::autoPtr<Foam::PackedList<nBits>>
Foam::PackedList<nBits>::clone() const
{
    return autoPtr<PackedList<nBits>>::New(*this);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<unsigned nBits>
inline void Foam::PackedList<nBits>::checkIndex(const label i) const
{
    if (!size_)
    {
        FatalErrorInFunction
            << "attempt to access element " << i << " from zero sized list"
            << abort(FatalError);
    }
    else if (i < 0 || i >= size_)
    {
        FatalErrorInFunction
            << "index " << i << " out of range [0," << size_ << ")"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Reference  * * * * * * * * * * * * * * * * //

template<unsigned nBits>
inline Foam::PackedList<nBits>::reference::reference
(
    PackedList* list,
    const label index
)
:
    ref_(list->StorageList::operator[](index / packing())),
    shift_(nBits * (index % packing()))
{}


template<unsigned nBits>
inline unsigned int Foam::PackedList<nBits>::reference::get() const
{
    return ((ref_ >> shift_) & max_value());
}


template<unsigned nBits>
inline bool Foam::PackedList<nBits>::reference::set(const unsigned int val)
{
    const unsigned int prev = ref_;

    const unsigned int mask = max_value() << shift_;

    if (val >= max_value())
    {
        // Overflow is max_value, fill everything
        ref_ |= mask;
    }
    else
    {
        ref_ &= ~mask;
        ref_ |= mask & (val << shift_);
    }

    return (prev != ref_);
}


template<unsigned nBits>
inline void Foam::PackedList<nBits>::reference::operator=
(
    const reference& other
)
{
    this->set(other.get());
}


template<unsigned nBits>
inline void Foam::PackedList<nBits>::reference::operator=
(
    const unsigned int val
)
{
    this->set(val);
}


template<unsigned nBits>
inline Foam::PackedList<nBits>::reference::operator unsigned int () const
{
    return this->get();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<unsigned nBits>
inline Foam::label Foam::PackedList<nBits>::size() const
{
    return size_;
}


template<unsigned nBits>
inline bool Foam::PackedList<nBits>::empty() const
{
    return !size_;
}


template<unsigned nBits>
inline void Foam::PackedList<nBits>::resize
(
    const label newSize,
    const unsigned int val
)
{
    reserve(newSize);

    const label oldSize = size_;
    size_ = newSize;

    if (size_ > oldSize)
    {
        // Fill new elements or newly exposed elements
        if (val)
        {
            // Fill value for complete segments
            unsigned int fill = val;

            if (val >= max_value())
            {
                // Fill everything
                fill = maskLower(packing());
            }
            else
            {
                for (unsigned int i = 1; i < packing(); ++i)
                {
                    fill |= (fill << nBits);
                }
            }

            // Fill in complete segments
            const label oldLen = packedLength(oldSize);
            const label newLen = packedLength(size_);
            for (label i=oldLen; i < newLen; ++i)
            {
                StorageList::operator[](i) = fill;
            }

            // Finish previous partial segment, preserve existing value
            {
                const unsigned int off = oldSize % packing();
                if (off)
                {
                    const unsigned int seg = oldSize / packing();
                    const unsigned int mask = maskLower(off);

                    StorageList::operator[](seg) &= mask;
                    StorageList::operator[](seg) |= ~mask & fill;
                }
            }

            // Mask off the (new) final partial segment
            {
                const unsigned int off = size_ % packing();
                if (off)
                {
                    const unsigned int seg = size_ / packing();

                    StorageList::operator[](seg) &= maskLower(off);
                }
            }
        }
    }
    else if (size_ < oldSize)
    {
        // Resize shrinking
        // - clear newly exposed elements

        // Fill in complete segments
        const label oldLen = packedLength(oldSize);
        const label newLen = packedLength(size_);
        for (label i=newLen; i < oldLen; ++i)
        {
            StorageList::operator[](i) = 0u;
        }

        // Mask off the final partial segment
        {
            const unsigned int off = size_ % packing();
            if (off)
            {
                const unsigned int seg = size_ / packing();

                StorageList::operator[](seg) &= maskLower(off);
            }
        }
    }
}


template<unsigned nBits>
inline void Foam::PackedList<nBits>::setSize
(
    const label newSize,
    const unsigned int val
)
{
    resize(newSize, val);
}


template<unsigned nBits>
inline Foam::label Foam::PackedList<nBits>::capacity() const
{
    return packing() * StorageList::size();
}


template<unsigned nBits>
inline void Foam::PackedList<nBits>::setCapacity(const label nElem)
{
    StorageList::setSize(packedLength(nElem), 0u);

    // Truncate addressed size too
    if (size_ > nElem)
    {
        size_ = nElem;

        // Mask off the final partial segment
        const unsigned int off = size_ % packing();
        if (off)
        {
            const unsigned int seg = size_ / packing();

            StorageList::operator[](seg) &= maskLower(off);
        }
    }
}


template<unsigned nBits>
inline void Foam::PackedList<nBits>::reserve(const label nElem)
{
    const label len = packedLength(nElem);

    // Allocate more capacity if necessary
    if (len > StorageList::size())
    {
        StorageList::setSize
        (
            max
            (
                len,
                // SizeInc=0, SizeMult=2, SizeDiv=1
                2 * StorageList::size()
            ),
            0u
        );
    }
}


template<unsigned nBits>
inline void Foam::PackedList<nBits>::reset()
{
    StorageList::operator=(0u);
}


template<unsigned nBits>
inline void Foam::PackedList<nBits>::clear()
{
    reset();
    size_ = 0;
}


template<unsigned nBits>
inline void Foam::PackedList<nBits>::clearStorage()
{
    StorageList::clear();
    size_ = 0;
}


template<unsigned nBits>
inline void Foam::PackedList<nBits>::shrink()
{
    // Any unneeded space allocated?
    const label len = packedLength();
    if (len < StorageList::size())
    {
        StorageList::setSize(len);
    }
}


template<unsigned nBits>
inline Foam::List<unsigned int>& Foam::PackedList<nBits>::storage()
{
    return static_cast<StorageList&>(*this);
}


template<unsigned nBits>
inline const Foam::List<unsigned int>& Foam::PackedList<nBits>::storage() const
{
    return static_cast<const StorageList&>(*this);
}


template<unsigned nBits>
inline Foam::label Foam::PackedList<nBits>::packedLength() const
{
    return packedLength(size_);
}


template<unsigned nBits>
inline std::streamsize Foam::PackedList<nBits>::byteSize() const
{
    return packedLength() * sizeof(StorageType);
}


template<unsigned nBits>
inline void Foam::PackedList<nBits>::transfer(PackedList<nBits>& lst)
{
    size_ = lst.size_;
    lst.size_ = 0;

    StorageList::transfer(lst);
}


template<unsigned nBits>
inline unsigned int Foam::PackedList<nBits>::get(const label i) const
{
    // Lazy evaluation - return 0 for out-of-range
    if (i < 0 || i >= size_)
    {
        return 0u;
    }

    return reference(const_cast<PackedList<nBits>*>(this), i).get();
}


template<unsigned nBits>
inline bool Foam::PackedList<nBits>::set
(
    const label i,
    const unsigned int val
)
{
    if (i < 0)
    {
        // Lazy evaluation - ignore out-of-bounds
        return false;
    }
    else if (i >= size_)
    {
        if (!val)
        {
            // Same as unset out-of-bounds = noop
            return false;
        }

        // Lazy evaluation - increase size on assignment
        resize(i + 1);
    }

    return reference(this, i).set(val);
}


template<unsigned nBits>
inline bool Foam::PackedList<nBits>::unset(const label i)
{
    // Unset out-of-bounds = noop
    if (i < 0 || i >= size_)
    {
        return false;
    }

    return reference(this, i).set(0u);
}


template<unsigned nBits>
inline Foam::PackedList<nBits>&
Foam::PackedList<nBits>::append(const unsigned int val)
{
    const label idx = size_;
    reserve(idx + 1);
    size_++;

    reference(this, idx).set(val);
    return *this;
}


template<unsigned nBits>
inline unsigned int Foam::PackedList<nBits>::remove()
{
    // Location of last element and simultaneously the new size
    const label idx = size_ - 1;

    if (idx < 0)
    {
        FatalErrorInFunction
            << "List is empty" << abort(FatalError);
    }

    const unsigned int val = reference(this, idx).get();
    resize(idx);

    return val;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<unsigned nBits>
inline unsigned int Foam::PackedList<nBits>::operator[](const label i) const
{
    return get(i);
}


template<unsigned nBits>
inline typename Foam::PackedList<nBits>::reference
Foam::PackedList<nBits>::operator[](const label i)
{
    // Leave enabled during testing period (MAR-2018)
    // #ifdef FULLDEBUG
    checkIndex(i);
    // #endif
    return reference(this, i);
}


template<unsigned nBits>
inline void Foam::PackedList<nBits>::operator=(const unsigned int val)
{
    const label packLen = packedLength();

    if (val && size_)
    {
        unsigned int fill = val;

        if (val >= max_value())
        {
            // Fill everything
            fill = maskLower(packing());
        }
        else
        {
            for (unsigned int i = 1; i < packing(); ++i)
            {
                fill |= (fill << nBits);
            }
        }

        for (label i=0; i < packLen; ++i)
        {
            StorageList::operator[](i) = fill;
        }

        // Mask off the final partial segment
        {
            const unsigned int off = size_ % packing();
            if (off)
            {
                const unsigned int seg = size_ / packing();

                StorageList::operator[](seg) &= maskLower(off);
            }
        }
    }
    else
    {
        for (label i=0; i < packLen; ++i)
        {
            StorageList::operator[](i) = 0u;
        }
    }
}


// ************************************************************************* //
