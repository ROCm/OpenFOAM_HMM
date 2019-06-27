/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include <algorithm>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline constexpr Foam::sliceRange::sliceRange() noexcept
:
    start_(0),
    size_(0),
    stride_(0)
{}


inline Foam::sliceRange::sliceRange
(
    const label start,
    const label size,
    const label stride
) noexcept
:
    start_(start),
    size_(std::max(label(0),size)),    // No negative size
    stride_(std::max(label(0),stride)) // No negative stride
{}


// * * * * * * * * * * * * * * * * Iterators * * * * * * * * * * * * * * * * //

inline Foam::sliceRange::indexer::indexer
(
    const sliceRange* range,
    const label i
)
:
    stride_(range->stride()),
    value_
    (
        range->start()
      + ((i < 0 || i > range->size()) ? range->size() : i) * stride_
    )
{}


inline void Foam::sliceRange::indexer::next() noexcept
{
    value_ += stride_;
}


inline void Foam::sliceRange::indexer::next(const label n) noexcept
{
    value_ += (n * stride_);
}


inline bool Foam::sliceRange::indexer::equals
(
    const indexer& other
) const noexcept
{
    return (value_ == other.value_);
}


inline Foam::label Foam::sliceRange::indexer::operator*() const noexcept
{
    return value_;
}


inline Foam::label Foam::sliceRange::indexer::operator()()
{
    const label old(value_);
    next();
    return old;
}


inline Foam::sliceRange::const_iterator::const_iterator
(
    const sliceRange* range,
    const label i
)
:
    indexer(range, i)
{}


inline Foam::sliceRange::const_iterator&
Foam::sliceRange::const_iterator::operator++() noexcept
{
    next();
    return *this;
}


inline Foam::sliceRange::const_iterator&
Foam::sliceRange::const_iterator::operator+=(const label n) noexcept
{
    next(n);
    return *this;
}


inline Foam::sliceRange::const_iterator&
Foam::sliceRange::const_iterator::operator--() noexcept
{
    next(-1);
    return *this;
}


inline Foam::sliceRange::const_iterator&
Foam::sliceRange::const_iterator::operator-=(const label n) noexcept
{
    next(-n);
    return *this;
}


inline bool Foam::sliceRange::const_iterator::operator==
(
    const const_iterator& iter
) const noexcept
{
    return equals(iter);
}


inline bool Foam::sliceRange::const_iterator::operator!=
(
    const const_iterator& iter
) const noexcept
{
    return !equals(iter);
}


inline Foam::sliceRange::const_iterator Foam::sliceRange::begin() const
{
    return const_iterator(this, 0);
}


inline Foam::sliceRange::const_iterator Foam::sliceRange::cbegin() const
{
    return const_iterator(this, 0);
}


inline const Foam::sliceRange::const_iterator Foam::sliceRange::end() const
{
    return const_iterator(this, -1);
}


inline const Foam::sliceRange::const_iterator Foam::sliceRange::cend() const
{
    return const_iterator(this, -1);
}


inline Foam::sliceRange::indexer Foam::sliceRange::generator() const
{
    return indexer(this, 0);
}


inline Foam::sliceRange::const_iterator
Foam::sliceRange::at(const label i) const
{
    return const_iterator(this, i);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline Foam::label Foam::sliceRange::operator[](const label i) const noexcept
{
    return start_ + stride_ * i;
}


// ************************************************************************* //
