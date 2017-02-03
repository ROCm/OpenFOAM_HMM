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
    const labelRange& range,
    const bool endIter
)
:
   range_(range),
   index_(endIter ? range_.size() : 0)
{}


inline bool Foam::labelRange::const_iterator::operator==
(
    const const_iterator& iter
) const
{
    return (this->index_ == iter.index_);
}


inline bool Foam::labelRange::const_iterator::operator!=
(
    const const_iterator& iter
) const
{
    return (this->index_ != iter.index_);
}


inline Foam::label Foam::labelRange::const_iterator::operator*() const
{
    return range_[index_];
}


inline Foam::label Foam::labelRange::const_iterator::operator()() const
{
    return range_[index_];
}


inline Foam::labelRange::const_iterator&
Foam::labelRange::const_iterator::operator++()
{
    ++index_;
    return *this;
}


inline Foam::labelRange::const_iterator
Foam::labelRange::const_iterator::operator++(int)
{
    const_iterator old = *this;
    ++index_;
    return old;
}


inline Foam::labelRange::const_iterator Foam::labelRange::cbegin() const
{
    return const_iterator(*this);
}


inline const Foam::labelRange::const_iterator Foam::labelRange::cend() const
{
    return const_iterator(*this, true);
}


inline Foam::labelRange::const_iterator Foam::labelRange::begin() const
{
    return const_iterator(*this);
}


inline const Foam::labelRange::const_iterator Foam::labelRange::end() const
{
    return const_iterator(*this, true);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::labelRange::resize(const label n)
{
    setSize(n);
}


inline void Foam::labelRange::setSize(const label n)
{
    size_ = n;

    if (size_ < 0)
    {
        size_ = 0;
    }
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


inline Foam::label Foam::labelRange::first() const
{
    return start_;
}


inline Foam::label Foam::labelRange::last() const
{
    return start_ + size_ - 1;
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


inline bool Foam::labelRange::contains(const label value) const
{
    return value >= this->first() && value <= this->last();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline Foam::label Foam::labelRange::operator[](const label i) const
{
    return start_ + i;
}


inline bool Foam::labelRange::operator<(const labelRange& rhs) const
{
    return
    (
        this->first() < rhs.first()
     || (this->first() == rhs.first() && this->size() < rhs.size())
    );
}


inline bool Foam::labelRange::operator==(const labelRange& rhs) const
{
    return start_ == rhs.start_ && size_ == rhs.size_;
}


inline bool Foam::labelRange::operator!=(const labelRange& rhs) const
{
    return !(operator==(rhs));
}


// ************************************************************************* //
