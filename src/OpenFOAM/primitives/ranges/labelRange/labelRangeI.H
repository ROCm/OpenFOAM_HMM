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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::labelRange::labelRange()
:
    start_(0),
    size_(0)
{}


inline Foam::labelRange::labelRange
(
    const label start,
    const label size,
    const bool adjustStart
)
:
    start_(start),
    size_(size)
{
    if (adjustStart)
    {
        // Disallow invalid indices and sizes
        adjust();
    }
    else if (size_ < 0)
    {
        // Disallow invalid sizes
        size_ = 0;
    }
}


// * * * * * * * * * * * * * * * * Iterators * * * * * * * * * * * * * * * * //

inline Foam::labelRange::const_iterator::const_iterator
(
    const labelRange* range,
    const label i
)
:
    value_
    (
        range->start()
      + ((i < 0 || i > range->size()) ? range->size() : i)
    )
{}


inline Foam::label Foam::labelRange::const_iterator::operator*() const
{
    return value_;
}


inline Foam::labelRange::const_iterator&
Foam::labelRange::const_iterator::operator++()
{
    ++value_;
    return *this;
}


inline Foam::labelRange::const_iterator
Foam::labelRange::const_iterator::operator++(int)
{
    const_iterator old = *this;
    ++value_;
    return old;
}


inline bool Foam::labelRange::const_iterator::operator==
(
    const const_iterator& iter
) const
{
    return (this->value_ == iter.value_);
}


inline bool Foam::labelRange::const_iterator::operator!=
(
    const const_iterator& iter
) const
{
    return (this->value_ != iter.value_);
}


inline Foam::labelRange::const_iterator Foam::labelRange::begin() const
{
    return const_iterator(this, 0);
}


inline Foam::labelRange::const_iterator Foam::labelRange::cbegin() const
{
    return const_iterator(this, 0);
}


inline const Foam::labelRange::const_iterator Foam::labelRange::end() const
{
    return const_iterator(this, -1);
}


inline const Foam::labelRange::const_iterator Foam::labelRange::cend() const
{
    return const_iterator(this, -1);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::labelRange Foam::labelRange::identity(const label len)
{
    return labelRange(0, len);
}


inline void Foam::labelRange::setStart(const label i)
{
    start_ = i;
}


inline void Foam::labelRange::setSize(const label n)
{
    size_ = n;
    if (size_ < 0) size_ = 0;
}


inline void Foam::labelRange::resize(const label n)
{
    setSize(n);
}


inline void Foam::labelRange::decrement()
{
    --size_;
    if (size_ < 0) size_ = 0;
}


inline void Foam::labelRange::increment()
{
    ++size_;
}


inline void Foam::labelRange::clear()
{
    start_ = size_ = 0;
}


inline bool Foam::labelRange::empty() const
{
    return !size_;
}


inline bool Foam::labelRange::valid() const
{
    return size_;
}


inline Foam::label Foam::labelRange::size() const
{
    return size_;
}


inline Foam::label Foam::labelRange::start() const
{
    return start_;
}


inline Foam::label Foam::labelRange::first() const
{
    return start_;
}


inline Foam::label Foam::labelRange::last() const
{
    return start_ + size_ - 1;
}


inline Foam::label Foam::labelRange::before() const
{
    return start_ - 1;
}


inline Foam::label Foam::labelRange::after() const
{
    return start_ + size_;
}


inline bool Foam::labelRange::reset
(
    const label start,
    const label size,
    const bool adjustStart
)
{
    start_ = start;
    size_  = size;

    if (adjustStart)
    {
        // Disallow invalid indices and sizes
        adjust();
    }
    else if (size_ < 0)
    {
        // Disallow invalid sizes
        size_ = 0;
    }

    return size_;
}


inline bool Foam::labelRange::found(const label value) const
{
    return (value >= this->first() && value <= this->last());
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline Foam::label Foam::labelRange::operator[](const label localIndex) const
{
    return start_ + localIndex;
}


inline Foam::labelRange::const_iterator
Foam::labelRange::operator()(const label localIndex) const
{
    return const_iterator(this, localIndex);
}


inline Foam::label Foam::labelRange::operator++()
{
    return ++size_;
}


inline Foam::label Foam::labelRange::operator++(int)
{
    return size_++;
}


inline Foam::label Foam::labelRange::operator--()
{
    decrement();
    return size_;
}


inline Foam::label Foam::labelRange::operator--(int)
{
    const label old = size_;
    decrement();
    return old;
}


// ************************************************************************* //
