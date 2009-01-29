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
#include <climits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<unsigned nBits>
inline unsigned int Foam::PackedList<nBits>::max_bits()
{
    return sizeof(StorageType)*CHAR_BIT - 1;
}


template<unsigned nBits>
inline unsigned int Foam::PackedList<nBits>::max_value()
{
    return (1 << nBits) - 1;
}


template<unsigned nBits>
inline unsigned int Foam::PackedList<nBits>::packing()
{
    return sizeof(StorageType)*CHAR_BIT / nBits;
}


template<unsigned nBits>
inline unsigned int Foam::PackedList<nBits>::maskLower(unsigned offset)
{
    return (1 << (nBits * offset)) - 1;
}


template<unsigned nBits>
inline Foam::label Foam::PackedList<nBits>::packedLength(const label nElem)
{
    return (nElem + packing() - 1) / packing();
}


template<unsigned nBits>
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

template<unsigned nBits>
inline Foam::PackedList<nBits>::PackedList()
:
    StorageList(),
    size_(0)
{}


template<unsigned nBits>
inline Foam::PackedList<nBits>::PackedList(const label size)
:
    StorageList(packedLength(size), 0u),
    size_(size)
{}


template<unsigned nBits>
inline Foam::PackedList<nBits>::PackedList(const PackedList<nBits>& lst)
:
    StorageList(lst),
    size_(lst.size())
{}


template<unsigned nBits>
inline Foam::PackedList<nBits>::PackedList(const Xfer<PackedList<nBits> >& lst)
{
    transfer(lst());
}


template<unsigned nBits>
inline Foam::autoPtr<Foam::PackedList<nBits> >
Foam::PackedList<nBits>::clone() const
{
    return autoPtr<PackedList<nBits> >(new PackedList<nBits>(*this));
}


// * * * * * * * * * * * * * * * * Iterators * * * * * * * * * * * * * * * * //

// iteratorBase

template<unsigned nBits>
inline Foam::PackedList<nBits>::iteratorBase::iteratorBase()
:
    list_(0),
    index_(0),
    offset_(0)
{}


template<unsigned nBits>
inline Foam::PackedList<nBits>::iteratorBase::iteratorBase
(
    const iteratorBase& iter
)
:
    list_(iter.list_),
    index_(iter.index_),
    offset_(iter.offset_)
{}


template<unsigned nBits>
inline Foam::PackedList<nBits>::iteratorBase::iteratorBase
(
    const PackedList<nBits>* lst,
    const label i
)
:
    list_(const_cast<PackedList<nBits>*>(lst)),
    index_(i / packing()),
    offset_(i % packing())
{}


template<unsigned nBits>
inline unsigned int
Foam::PackedList<nBits>::iteratorBase::get() const
{
    const unsigned int& stored = list_->StorageList::operator[](index_);
    return (stored >> (nBits * offset_)) & max_value();
}


template<unsigned nBits>
inline bool
Foam::PackedList<nBits>::iteratorBase::set(const unsigned int val)
{
    unsigned int& stored = list_->StorageList::operator[](index_);
    const unsigned int prev = stored;

    const unsigned int startBit = nBits * offset_;
    const unsigned int maskNew  = max_value() << startBit;

    if (val & ~max_value())
    {
#       ifdef DEBUGList
        FatalErrorIn("PackedList<T>::iteratorBase::set(const unsigned int)")
            << "value " << label(val)
w            << " out-of-range 0 ... " << label(max_value())
            << " representable by " << nBits << " bits"
            << abort(FatalError);
#       endif

        // treat overflow as max_value
        stored |= maskNew;
    }
    else
    {
        stored = (stored & ~maskNew) | (maskNew & (val << startBit));
    }

    return prev != stored;
}



template<unsigned nBits>
inline void
Foam::PackedList<nBits>::iteratorBase::incr()
{
    offset_++;
    if (offset_ >= packing())
    {
        offset_ = 0;
        index_++;
    }
}


template<unsigned nBits>
inline void
Foam::PackedList<nBits>::iteratorBase::decr()
{
    if (!offset_)
    {
        offset_ = packing();
        index_--;
    }
    offset_--;
}


template<unsigned nBits>
inline void
Foam::PackedList<nBits>::iteratorBase::seek
(
    const iteratorBase& iter
)
{
    list_   = iter.list_;
    index_  = iter.index_;
    offset_ = iter.offset_;

    this->validate();
}


template<unsigned nBits>
inline bool
Foam::PackedList<nBits>::iteratorBase::valid() const
{
    label elemI = offset_ + index_ * packing();
    return (elemI < list_->size_);
}


template<unsigned nBits>
inline bool
Foam::PackedList<nBits>::iteratorBase::validate()
{
    // avoid going past end()
    unsigned endIdx = list_->size_ / packing();
    unsigned endOff = list_->size_ % packing();

    if (index_ > endIdx || (index_ == endIdx && offset_ > endOff))
    {
        index_  = endIdx;
        offset_ = endOff;

        return false;
    }
    else
    {
        return true;
    }
}


template<unsigned nBits>
inline bool Foam::PackedList<nBits>::iteratorBase::operator==
(
    const iteratorBase& iter
) const
{
    return this->index_ == iter.index_ && this->offset_ == iter.offset_;
}


template<unsigned nBits>
inline bool Foam::PackedList<nBits>::iteratorBase::operator!=
(
    const iteratorBase& iter
) const
{
    return this->index_ != iter.index_ || this->offset_ != iter.offset_;
}


template<unsigned nBits>
inline unsigned int
Foam::PackedList<nBits>::iteratorBase::operator=(const iteratorBase& iter)
{
    const unsigned int val = iter.get();
    this->set(val);
    return val;
}


template<unsigned nBits>
inline unsigned int
Foam::PackedList<nBits>::iteratorBase::operator=(const unsigned int val)
{
#   ifdef DEBUGList
    // lazy evaluation would be nice to keep, but really slows things down
    label minsize = 1 + offset_ + index_ * packing();
    if (minsize > list_->size_)
    {
        list_->resize(minsize);
    }
#endif

    this->set(val);
    return val;
}


template<unsigned nBits>
inline Foam::PackedList<nBits>::iteratorBase::operator
unsigned int () const
{
#   ifdef DEBUGList
    // lazy evaluation would be nice to keep, but really slows things down
    label minsize = 1 + offset_ + index_ * packing();
    if (minsize > list_->size_)
    {
        return 0;
    }
#endif

    return this->get();
}


// const_iterator, iterator

template<unsigned nBits>
inline Foam::PackedList<nBits>::iterator::iterator()
:
    iteratorBase()
{}


template<unsigned nBits>
inline Foam::PackedList<nBits>::const_iterator::const_iterator()
:
    iteratorBase()
{}


template<unsigned nBits>
inline Foam::PackedList<nBits>::iterator::iterator
(
    const iteratorBase& iter
)
:
    iteratorBase(iter)
{
    this->validate();
}


template<unsigned nBits>
inline Foam::PackedList<nBits>::const_iterator::const_iterator
(
    const iteratorBase& iter
)
:
    iteratorBase(iter)
{
    this->validate();
}


template<unsigned nBits>
inline Foam::PackedList<nBits>::iterator::iterator
(
    const PackedList<nBits>* lst,
    const label i
)
:
    iteratorBase(lst, i)
{}


template<unsigned nBits>
inline Foam::PackedList<nBits>::const_iterator::const_iterator
(
    const PackedList<nBits>* lst,
    const label i
)
:
    iteratorBase(lst, i)
{}


template<unsigned nBits>
inline Foam::PackedList<nBits>::const_iterator::const_iterator
(
    const iterator& iter
)
:
    iteratorBase(static_cast<const iteratorBase&>(iter))
{}


template<unsigned nBits>
inline typename Foam::PackedList<nBits>::iterator&
Foam::PackedList<nBits>::iterator::operator=(const iteratorBase& iter)
{
    this->seek(iter);
    return *this;
}


template<unsigned nBits>
inline typename Foam::PackedList<nBits>::const_iterator&
Foam::PackedList<nBits>::const_iterator::operator=(const iteratorBase& iter)
{
    this->seek(iter);
    return *this;
}


template<unsigned nBits>
inline typename Foam::PackedList<nBits>::iterator&
Foam::PackedList<nBits>::iterator::operator++()
{
    this->incr();
    return *this;
}


template<unsigned nBits>
inline typename Foam::PackedList<nBits>::const_iterator&
Foam::PackedList<nBits>::const_iterator::operator++()
{
    this->incr();
    return *this;
}


template<unsigned nBits>
inline typename Foam::PackedList<nBits>::iterator
Foam::PackedList<nBits>::iterator::operator++(int)
{
    iterator old = *this;
    this->incr();
    return old;
}


template<unsigned nBits>
inline typename Foam::PackedList<nBits>::const_iterator
Foam::PackedList<nBits>::const_iterator::operator++(int)
{
    const_iterator old = *this;
    this->incr();
    return old;
}


template<unsigned nBits>
inline typename Foam::PackedList<nBits>::iterator&
Foam::PackedList<nBits>::iterator::operator--()
{
    this->decr();
    return *this;
}


template<unsigned nBits>
inline typename Foam::PackedList<nBits>::const_iterator&
Foam::PackedList<nBits>::const_iterator::operator--()
{
    this->decr();
    return *this;
}


template<unsigned nBits>
inline typename Foam::PackedList<nBits>::iterator
Foam::PackedList<nBits>::iterator::operator--(int)
{
    iterator old = *this;
    this->decr();
    return old;
}


template<unsigned nBits>
inline typename Foam::PackedList<nBits>::const_iterator
Foam::PackedList<nBits>::const_iterator::operator--(int)
{
    const_iterator old = *this;
    this->decr();
    return old;
}


template<unsigned nBits>
inline typename Foam::PackedList<nBits>::iteratorBase&
Foam::PackedList<nBits>::iterator::operator*()
{
    return static_cast<iteratorBase&>(*this);
}


template<unsigned nBits>
inline typename Foam::PackedList<nBits>::iteratorBase&
Foam::PackedList<nBits>::iterator::operator()()
{
    return static_cast<iteratorBase&>(*this);
}


template<unsigned nBits>
inline unsigned int
Foam::PackedList<nBits>::const_iterator::operator*() const
{
    return this->get();
}


template<unsigned nBits>
inline unsigned int
Foam::PackedList<nBits>::const_iterator::operator()() const
{
    return this->get();
}


template<unsigned nBits>
inline typename Foam::PackedList<nBits>::iterator
Foam::PackedList<nBits>::begin()
{
    return iterator(this, 0);
}


template<unsigned nBits>
inline typename Foam::PackedList<nBits>::const_iterator
Foam::PackedList<nBits>::begin() const
{
    return const_iterator(this, 0);
}


template<unsigned nBits>
inline typename Foam::PackedList<nBits>::const_iterator
Foam::PackedList<nBits>::cbegin() const
{
    return const_iterator(this, 0);
}


template<unsigned nBits>
inline typename Foam::PackedList<nBits>::iterator
Foam::PackedList<nBits>::end()
{
    return iterator(this, this->size());
}


template<unsigned nBits>
inline typename Foam::PackedList<nBits>::const_iterator
Foam::PackedList<nBits>::end() const
{
    return const_iterator(this, this->size());
}


template<unsigned nBits>
inline typename Foam::PackedList<nBits>::const_iterator
Foam::PackedList<nBits>::cend() const
{
    return const_iterator(this, this->size());
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
    const label nElem,
    const unsigned int& val
)
{
    reserve(nElem);

    if (nElem > size_)
    {
        // fill new elements or newly exposed elements
        if (size_)
        {
            // fill value for complete chunks
            unsigned int fill = val;

            if (fill & ~max_value())
            {
#               ifdef DEBUGList
                FatalErrorIn("PackedList<T>::resize(label, const unsigned int)")
                    << "value " << label(val)
                    << " out-of-range 0 ... " << label(max_value())
                    << " representable by " << nBits << " bits"
                    << abort(FatalError);
#               endif

                // treat overflow as max_value, fill everything
                fill = ~0;
            }
            else
            {
                for (unsigned int i = 1; i < packing(); ++i)
                {
                    fill |= (fill << nBits);
                }
            }

            unsigned begIdx = size_ / packing();
            unsigned begOff = size_ % packing();
            unsigned endIdx = nElem / packing();

            // partial chunk, preserve existing value
            if (begOff)
            {
                unsigned int maskOld = maskLower(begOff);

                StorageList::operator[](begIdx) &= maskOld;
                StorageList::operator[](begIdx) |= ~maskOld & fill;

                // continue with the next chunk
                begIdx++;
            }

            // fill in complete elements
            while (begIdx < endIdx)
            {
                StorageList::operator[](begIdx++) = fill;
            }
        }
        else
        {
            // no original size - simply flood-fill
            operator=(val);
        }
    }

    size_ = nElem;
}


template<unsigned nBits>
inline void Foam::PackedList<nBits>::setSize
(
    const label newSize,
    const unsigned int& val
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

    // truncate addressed size too?
    if (size_ > nElem)
    {
        size_ = nElem;
    }
}


template<unsigned nBits>
inline void Foam::PackedList<nBits>::reserve
(
    const label nElem
)
{
    label len = packedLength(nElem);

    // need more capacity?
    if (len > StorageList::size())
    {
        StorageList::setSize(len, 0u);
    }
}


template<unsigned nBits>
inline void Foam::PackedList<nBits>::clear()
{
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
    label len = packedLength(size_);

    // we have unused space?
    if (len < StorageList::size())
    {
        StorageList::setSize(len);
    }
}

template<unsigned nBits>
inline Foam::List<unsigned int>&
Foam::PackedList<nBits>::storage()
{
    return static_cast<StorageList&>(*this);
}


template<unsigned nBits>
inline const Foam::List<unsigned int>&
Foam::PackedList<nBits>::storage() const
{
    return static_cast<const StorageList&>(*this);
}


template<unsigned nBits>
inline void Foam::PackedList<nBits>::transfer(PackedList<nBits>& lst)
{
    size_ = lst.size_;
    lst.size_ = 0;

    StorageList::transfer(lst);
}


template<unsigned nBits>
inline Foam::Xfer< Foam::PackedList<nBits> >
Foam::PackedList<nBits>::xfer()
{
    return xferMove(*this);
}


template<unsigned nBits>
inline unsigned int Foam::PackedList<nBits>::get(const label i) const
{
#   ifdef DEBUGList
    checkIndex(i);
#   endif

    return iteratorBase(this, i).get();
}


template<unsigned nBits>
inline unsigned int Foam::PackedList<nBits>::operator[](const label i) const
{
    if (i < size_)
    {
        return iteratorBase(this, i).get();
    }
    else
    {
        return 0;
    }
}


template<unsigned nBits>
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
    return iteratorBase(this, i).set(val);
}


template<unsigned nBits>
inline void Foam::PackedList<nBits>::append(const unsigned int val)
{
    label elemI = size_;
    reserve(elemI + 1);
    size_++;

    iteratorBase(this, elemI).set(val);
}


template<unsigned nBits>
inline unsigned int Foam::PackedList<nBits>::remove()
{
    if (!size_)
    {
        FatalErrorIn
        (
            "Foam::PackedList<nBits>::remove()"
        )   << "List is empty" << abort(FatalError);
    }

    label elemI = size_ - 1;
    const unsigned int val = iteratorBase(this, elemI).get();
    resize(elemI);

    return val;
}


template<unsigned nBits>
inline typename Foam::PackedList<nBits>::iteratorBase
Foam::PackedList<nBits>::operator[](const label i)
{
    return iteratorBase(this, i);
}


// specialization for nBits=1 isn't worth the bother
template<unsigned nBits>
inline void Foam::PackedList<nBits>::operator=(const unsigned int val)
{
    if (val)
    {
        unsigned int fill = val;

        if (fill & ~max_value())
        {
#           ifdef DEBUGList
            FatalErrorIn("PackedList<T>::operator=(const unsigned int)")
                << "value " << label(val)
                << " out-of-range 0 ... " << label(max_value())
                << " representable by " << nBits << " bits"
                << abort(FatalError);
#           endif

            // treat overflow as max_value
            fill = ~0;
        }
        else
        {
            for (unsigned int i = 1; i < packing(); ++i)
            {
                fill |= (fill << nBits);
            }
        }

        StorageList::operator=(fill);
    }
    else
    {
        StorageList::operator=(0u);
    }
}


// ************************************************************************* //
