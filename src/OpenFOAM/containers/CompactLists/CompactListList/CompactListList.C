/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "CompactListList.H"
#include "labelRange.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class T>
void Foam::CompactListList<T>::reportOverflowAndExit
(
    const label idx,
    const labelUList& localLens
)
{
    FatalErrorInFunction
        << "Overflow : sum of sizes exceeds labelMax ("
        << labelMax << ") after index " << idx;

    if (!localLens.empty())
    {
        FatalError << " of " << flatOutput(localLens);
    }

    FatalError
        << nl
        << "Please recompile with larger datatype for label." << nl
        << exit(FatalError);
}


template<class T>
template<class ListListType>
Foam::CompactListList<T> Foam::CompactListList<T>::packImpl
(
    const ListListType& lists,
    const bool checkOverflow
)
{
    CompactListList<T> compact;

    auto& newOffsets = compact.offsets_;
    auto& newValues = compact.values_;

    label total = 0;
    const label len = lists.size();

    if (len)
    {
        newOffsets.resize(len+1, Zero);

        for (label i = 0; i < len; ++i)
        {
            newOffsets[i] = total;
            total += lists[i].size();

            if (checkOverflow && total < newOffsets[i])
            {
                reportOverflowAndExit(i);
            }
        }
        newOffsets[len] = total;
    }

    if (total)
    {
        // Copy in the data
        newValues.resize(total);

        auto outIter = newValues.begin();

        for (const auto& list : lists)
        {
            forAll(list, i)
            {
                *outIter = list[i];
                ++outIter;
            }
        }
    }

    return compact;
}


template<class T>
template<class SubListType>
Foam::CompactListList<T> Foam::CompactListList<T>::pack
(
    const UList<SubListType>& lists,
    const bool checkOverflow
)
{
    return CompactListList<T>::packImpl<UList<SubListType>>
    (
        lists,
        checkOverflow
    );
}


template<class T>
template<class SubListType, class Addr>
Foam::CompactListList<T> Foam::CompactListList<T>::pack
(
    const IndirectListBase<SubListType, Addr>& lists,
    const bool checkOverflow
)
{
    return CompactListList<T>::packImpl<IndirectListBase<SubListType, Addr>>
    (
        lists,
        checkOverflow
    );
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
Foam::CompactListList<T>::CompactListList(const labelUList& listSizes)
{
    const label len = listSizes.size();

    if (len)
    {
        offsets_.resize(len+1);

        label total = 0;
        for (label i = 0; i < len; ++i)
        {
            offsets_[i] = total;
            total += listSizes[i];

#ifdef FULLDEBUG
            if (total < offsets_[i])
            {
                reportOverflowAndExit(i, listSizes);
            }
#endif
        }

        offsets_[len] = total;
        values_.resize(total);
    }
}


template<class T>
Foam::CompactListList<T>::CompactListList
(
    const labelUList& listSizes,
    const T& val
)
{
    const label len = listSizes.size();

    if (len)
    {
        offsets_.resize(len+1);

        label total = 0;
        for (label i = 0; i < len; ++i)
        {
            offsets_[i] = total;
            total += listSizes[i];

#ifdef FULLDEBUG
            if (total < offsets_[i])
            {
                reportOverflowAndExit(i, listSizes);
            }
#endif
        }

        offsets_[len] = total;
        values_.resize(total, val);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
Foam::label Foam::CompactListList<T>::maxNonLocalSize(const label rowi) const
{
    const label len = (offsets_.size() - 1);

    if (len < 1)
    {
        return 0;
    }

    label maxLen = 0;

    for (label i=0; i < len; ++i)
    {
        if (i != rowi)
        {
            const label localLen = (offsets_[i+1] - offsets_[i]);
            maxLen = max(maxLen, localLen);
        }
    }

    return maxLen;
}


template<class T>
std::streamsize Foam::CompactListList<T>::byteSize() const
{
    if (!is_contiguous<T>::value)
    {
        FatalErrorInFunction
            << "Invalid for non-contiguous data types"
            << abort(FatalError);
    }
    return this->size_bytes();
}


template<class T>
Foam::labelRange Foam::CompactListList<T>::range(const label i) const
{
    return labelRange(offsets_[i], offsets_[i+1] - offsets_[i]);
}


template<class T>
Foam::List<Foam::labelRange>
Foam::CompactListList<T>::ranges() const
{
    List<labelRange> values;

    const label len = (offsets_.size() - 1);

    if (len < 1)
    {
        return values;
    }

    values.resize(len);

    for (label i=0; i < len; ++i)
    {
        values[i].reset(offsets_[i], (offsets_[i+1] - offsets_[i]));
    }

    return values;
}


template<class T>
void Foam::CompactListList<T>::resize(const labelUList& listSizes)
{
    const label len = listSizes.size();

    if (len)
    {
        offsets_.resize(len+1);

        label total = 0;
        for (label i = 0; i < len; ++i)
        {
            offsets_[i] = total;
            total += listSizes[i];
#if 0
            if (checkOverflow && total < offsets_[i])
            {
                reportOverflowAndExit(i, listSizes);
            }
#endif
        }

        offsets_[len] = total;
        values_.resize(total);
    }
    else
    {
        clear();
    }
}


template<class T>
void Foam::CompactListList<T>::setLocalSize(const label rowi, const label len)
{
    if (rowi >= 0 && rowi+1 < offsets_.size() && len >= 0)
    {
        const label delta = (len - (offsets_[rowi+1] - offsets_[rowi]));

        // TBD: additional overflow check
        if (delta)
        {
            for (label i = rowi+1; i < offsets_.size(); ++i)
            {
                offsets_[i] += delta;
            }
        }
    }
}


template<class T>
Foam::labelList Foam::CompactListList<T>::localSizes() const
{
    labelList values;

    const label len = (offsets_.size() - 1);

    if (len < 1)
    {
        return values;
    }

    values.resize(len);

    for (label i=0; i < len; ++i)
    {
        values[i] = offsets_[i+1] - offsets_[i];
    }

    return values;
}


template<class T>
void Foam::CompactListList<T>::swap
(
    CompactListList<T>& other
)
{
    if (this == &other)
    {
        return;  // Self-swap is a no-op
    }

    offsets_.swap(other.offsets_);
    values_.swap(other.values_);
}


template<class T>
void Foam::CompactListList<T>::transfer
(
    CompactListList<T>& list
)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    offsets_.transfer(list.offsets_);
    values_.transfer(list.values_);
}


template<class T>
template<class SubListType>
Foam::List<SubListType>
Foam::CompactListList<T>::unpack() const
{
    List<SubListType> lists(size());

    forAll(lists, i)
    {
        lists[i] = SubListType(this->localList(i));
    }

    return lists;
}


template<class T>
template<class SubListType>
Foam::List<SubListType>
Foam::CompactListList<T>::unpack(const labelRange& range) const
{
    List<SubListType> lists(range.size());

    auto outIter = lists.begin();

    for (const label i : range)
    {
        *outIter = SubListType(this->localList(i));
        ++outIter;
    }

    return lists;
}


// ************************************************************************* //
