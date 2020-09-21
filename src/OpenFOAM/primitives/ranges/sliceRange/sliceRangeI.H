/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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


inline constexpr Foam::sliceRange::sliceRange
(
    const label beg,
    const label len,
    const label stride
) noexcept
:
    start_(beg),
    size_(len),
    stride_(stride)
{}


// * * * * * * * * * * * * * * * * Iterators * * * * * * * * * * * * * * * * //

inline constexpr Foam::sliceRange::indexer::indexer() noexcept
:
    stride_(1),
    value_(0)
{}


inline constexpr Foam::sliceRange::indexer::indexer
(
    const label val,
    const label stride
) noexcept
:
    stride_(stride),
    value_(val)
{}


inline Foam::label
Foam::sliceRange::indexer::
operator()()
{
    const label old(value_);
    next();
    return old;
}


inline Foam::sliceRange::const_iterator&
Foam::sliceRange::const_iterator::
operator++() noexcept
{
    next();
    return *this;
}


inline Foam::sliceRange::const_iterator&
Foam::sliceRange::const_iterator::
operator+=(const label n) noexcept
{
    next(n);
    return *this;
}


inline Foam::sliceRange::const_iterator&
Foam::sliceRange::const_iterator::
operator--() noexcept
{
    prev();
    return *this;
}


inline Foam::sliceRange::const_iterator&
Foam::sliceRange::const_iterator::
operator-=(const label n) noexcept
{
    prev(n);
    return *this;
}


inline Foam::sliceRange::const_reverse_iterator&
Foam::sliceRange::const_reverse_iterator::
operator++() noexcept
{
    prev();
    return *this;
}


inline Foam::sliceRange::const_reverse_iterator&
Foam::sliceRange::const_reverse_iterator::
operator+=(const label n) noexcept
{
    prev(n);
    return *this;
}


inline Foam::sliceRange::const_reverse_iterator&
Foam::sliceRange::const_reverse_iterator::
operator--() noexcept
{
    next();
    return *this;
}


inline Foam::sliceRange::const_reverse_iterator&
Foam::sliceRange::const_reverse_iterator::
operator-=(const label n) noexcept
{
    next(n);
    return *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::sliceRange::indexer
Foam::sliceRange::generator() const
{
    return indexer(start_, stride_);
}


inline Foam::sliceRange::const_iterator
Foam::sliceRange::at(const label i) const
{
    return
        const_iterator
        (
            start_ + ((i < 0 || i > size_) ? size_ : i) * stride_,
            stride_
        );
}


inline Foam::sliceRange::const_iterator
Foam::sliceRange::begin() const
{
    return const_iterator(start_, stride_);
}


inline Foam::sliceRange::const_iterator
Foam::sliceRange::cbegin() const
{
    return const_iterator(start_, stride_);
}


inline Foam::sliceRange::const_iterator
Foam::sliceRange::end() const
{
    return const_iterator(start_ + size_*stride_, stride_);
}


inline Foam::sliceRange::const_iterator
Foam::sliceRange::cend() const
{
    return const_iterator(start_ + size_*stride_, stride_);
}


inline Foam::sliceRange::const_reverse_iterator
Foam::sliceRange::rbegin() const
{
    return const_reverse_iterator(start_ + (size_-1)*stride_, stride_);
}


inline Foam::sliceRange::const_reverse_iterator
Foam::sliceRange::crbegin() const
{
    return const_reverse_iterator(start_ + (size_-1)*stride_, stride_);
}


inline Foam::sliceRange::const_reverse_iterator
Foam::sliceRange::rend() const
{
    return const_reverse_iterator(start_ - stride_, stride_);
}


inline Foam::sliceRange::const_reverse_iterator
Foam::sliceRange::crend() const
{
    return const_reverse_iterator(start_ - stride_, stride_);
}


// ************************************************************************* //
