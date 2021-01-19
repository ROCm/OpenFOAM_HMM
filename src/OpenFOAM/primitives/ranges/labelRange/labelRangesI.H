/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

inline Foam::labelRanges::labelRanges(const label nElem)
:
    DynamicList<labelRange>(nElem)
{}


// * * * * * * * * * * * * * * * * Iterators * * * * * * * * * * * * * * * * //

inline constexpr Foam::labelRanges::const_iterator::
const_iterator
(
    const labelRanges* list,
    const label i
) noexcept
:
    list_(list),
    index_(i),
    subIndex_(0)
{}


inline Foam::label Foam::labelRanges::const_iterator::
operator*() const
{
    return (*list_)[index_][subIndex_];
}


inline Foam::labelRanges::const_iterator&
Foam::labelRanges::const_iterator::
operator++()
{
    if (++subIndex_ >= (*list_)[index_].size())
    {
        // Move to the next range
        ++index_;
        subIndex_ = 0;
    }

    return *this;
}


inline Foam::labelRanges::const_iterator
Foam::labelRanges::const_iterator::
operator++(int)
{
    const_iterator old(*this);
    this->operator++();
    return old;
}


inline constexpr bool
Foam::labelRanges::const_iterator::
operator==
(
    const const_iterator& iter
) const noexcept
{
    return
    (
        index_    == iter.index_
     && subIndex_ == iter.subIndex_
    );
}


inline constexpr bool
Foam::labelRanges::const_iterator::
operator!=
(
    const const_iterator& iter
) const noexcept
{
    return !(*this == iter);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::labelRanges::const_iterator
Foam::labelRanges::cbegin() const noexcept
{
    return const_iterator(this);
}


inline const Foam::labelRanges::const_iterator
Foam::labelRanges::cend() const noexcept
{
    return const_iterator(this, this->size());
}


inline Foam::labelRanges::const_iterator
Foam::labelRanges::begin() const noexcept
{
    return const_iterator(this);
}


inline const Foam::labelRanges::const_iterator
Foam::labelRanges::end() const noexcept
{
    return const_iterator(this, this->size());
}


inline bool Foam::labelRanges::found(const label value) const
{
    for
    (
        const labelRange& range
      : static_cast<const StorageContainer&>(*this)
    )
    {
        if (range.found(value))
        {
            return true;
        }
    }

    return false;
}


// ************************************************************************* //
