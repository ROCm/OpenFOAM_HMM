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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline constexpr Foam::labelRange::labelRange() noexcept
:
    start_(0),
    size_(0)
{}


inline Foam::labelRange::labelRange
(
    const label start,
    const label size
) noexcept
:
    start_(start),
    size_(size)
{
    if (size_ < 0) size_ = 0;  // No negative sizes
}


inline Foam::labelRange::labelRange
(
    const label start,
    const label size,
    const bool adjustStart
) noexcept
:
    start_(start),
    size_(size)
{
    if (adjustStart)
    {
        adjust(); // Eliminate negative start, adjust size accordingly
    }
    else if (size_ < 0)
    {
        size_ = 0;  // No negative sizes
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


inline Foam::label Foam::labelRange::const_iterator::operator*() const noexcept
{
    return value_;
}


inline Foam::labelRange::const_iterator&
Foam::labelRange::const_iterator::operator++() noexcept
{
    ++value_;
    return *this;
}


inline Foam::labelRange::const_iterator
Foam::labelRange::const_iterator::operator++(int)
{
    const_iterator old(*this);
    ++value_;
    return old;
}


inline Foam::labelRange::const_iterator&
Foam::labelRange::const_iterator::operator+=(const label n) noexcept
{
    value_ += n;
    return *this;
}


inline Foam::labelRange::const_iterator&
Foam::labelRange::const_iterator::operator-=(const label n) noexcept
{
    value_ -= n;
    return *this;
}


inline bool Foam::labelRange::const_iterator::operator==
(
    const const_iterator& iter
) const noexcept
{
    return (value_ == iter.value_);
}


inline bool Foam::labelRange::const_iterator::operator!=
(
    const const_iterator& iter
) const noexcept
{
    return (value_ != iter.value_);
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


inline Foam::labelRange::const_iterator
Foam::labelRange::at(const label i) const
{
    return const_iterator(this, i);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::labelRange Foam::labelRange::identity
(
    const label len,
    const label start
) noexcept
{
    return labelRange(start, len);
}


inline void Foam::labelRange::setStart(const label i) noexcept
{
    start_ = i;
}


inline void Foam::labelRange::resize(const label n) noexcept
{
    size_ = n;
    if (size_ < 0) size_ = 0;  // No negative sizes
}


inline void Foam::labelRange::setSize(const label n) noexcept
{
    size_ = n;
    if (size_ < 0) size_ = 0;  // No negative sizes
}


inline void Foam::labelRange::clear() noexcept
{
    start_ = size_ = 0;
}


inline bool Foam::labelRange::empty() const noexcept
{
    return !size_;
}


inline bool Foam::labelRange::valid() const noexcept
{
    return size_;
}


inline Foam::label Foam::labelRange::size() const noexcept
{
    return size_;
}


inline Foam::label Foam::labelRange::start() const noexcept
{
    return start_;
}


inline Foam::label Foam::labelRange::before() const noexcept
{
    return start_ - 1;
}


inline Foam::label Foam::labelRange::first() const noexcept
{
    return start_;
}


inline Foam::label Foam::labelRange::min() const noexcept
{
    return start_;
}


inline Foam::label Foam::labelRange::last() const noexcept
{
    return start_ + size_ - 1;
}


inline Foam::label Foam::labelRange::max() const noexcept
{
    return start_ + size_ - 1;
}


inline Foam::label Foam::labelRange::after() const noexcept
{
    return start_ + size_;
}


inline bool Foam::labelRange::reset
(
    const label start,
    const label size
) noexcept
{
    start_ = start;
    size_  = size;

    if (size_ < 0) size_ = 0;  // No negative sizes

    return size_;
}


inline bool Foam::labelRange::reset
(
    const label start,
    const label size,
    const bool adjustStart
) noexcept
{
    start_ = start;
    size_  = size;

    if (adjustStart)
    {
        adjust(); // Eliminate negative start, adjust size accordingly
    }
    else if (size_ < 0)
    {
        size_ = 0;  // No negative sizes
    }

    return size_;
}


inline bool Foam::labelRange::found(const label value) const noexcept
{
    return (size_ && start_ <= value && (value - start_) < size_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline Foam::label Foam::labelRange::operator[](const label i) const noexcept
{
    return start_ + i;
}


inline bool Foam::labelRange::operator()(const label value) const noexcept
{
    return found(value);
}


inline Foam::label Foam::labelRange::operator++() noexcept
{
    return ++size_;
}


inline Foam::label Foam::labelRange::operator++(int) noexcept
{
    return size_++;
}


inline Foam::label Foam::labelRange::operator--() noexcept
{
    --size_;
    if (size_ < 0) size_ = 0;  // No negative sizes
    return size_;
}


inline Foam::label Foam::labelRange::operator--(int) noexcept
{
    const label old(size_);
    --size_;
    if (size_ < 0) size_ = 0;  // No negative sizes
    return old;
}


inline Foam::label Foam::labelRange::operator+=(const label n) noexcept
{
    size_ += n;
    return size_;
}


inline Foam::label Foam::labelRange::operator-=(const label n) noexcept
{
    size_ -= n;
    if (size_ < 0) size_ = 0;  // No negative sizes
    return size_;
}



// ************************************************************************* //
