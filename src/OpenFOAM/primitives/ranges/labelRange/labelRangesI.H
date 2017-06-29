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

inline Foam::labelRanges::labelRanges()
:
    StorageContainer()
{}


inline Foam::labelRanges::labelRanges(const label nElem)
:
    StorageContainer(nElem)
{}


// * * * * * * * * * * * * * * * * Iterators * * * * * * * * * * * * * * * * //

inline Foam::labelRanges::const_iterator::const_iterator
(
    const labelRanges* lst,
    const label i
)
:
    list_(lst),
    index_(i < 0 ? lst->size() : i),
    subIndex_(0)
{}


inline bool Foam::labelRanges::const_iterator::operator==
(
    const const_iterator& iter
) const
{
    return
    (
        this->index_    == iter.index_
     && this->subIndex_ == iter.subIndex_
    );
}


inline bool Foam::labelRanges::const_iterator::operator!=
(
    const const_iterator& iter
) const
{
    return !(this->operator==(iter));
}


inline Foam::label Foam::labelRanges::const_iterator::operator*()
{
    return list_->operator[](index_)[subIndex_];
}


inline Foam::labelRanges::const_iterator&
Foam::labelRanges::const_iterator::operator++()
{
    if (++subIndex_ >= list_->operator[](index_).size())
    {
        // Next sub-list
        ++index_;
        subIndex_ = 0;
    }

    return *this;
}


inline Foam::labelRanges::const_iterator
Foam::labelRanges::const_iterator::operator++(int)
{
    const_iterator old = *this;
    this->operator++();
    return old;
}


inline Foam::labelRanges::const_iterator Foam::labelRanges::cbegin() const
{
    return const_iterator(this, 0);
}


inline const Foam::labelRanges::const_iterator Foam::labelRanges::cend() const
{
    return const_iterator(this, -1);
}


inline Foam::labelRanges::const_iterator Foam::labelRanges::begin() const
{
    return const_iterator(this, 0);
}


inline const Foam::labelRanges::const_iterator Foam::labelRanges::end() const
{
    return const_iterator(this, -1);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::labelRanges::found(const label value) const
{
    forAll(*this, i)
    {
        if (StorageContainer::operator[](i).found(value))
        {
            return true;
        }
    }

    return false;
}


// ************************************************************************* //
