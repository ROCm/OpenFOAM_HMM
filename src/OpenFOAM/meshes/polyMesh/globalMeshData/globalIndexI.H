/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "ListOps.H"
#include "labelRange.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::globalIndex::globalIndex(const label localSize)
:
    globalIndex()
{
    reset(localSize);
}


inline Foam::globalIndex::globalIndex
(
    const label localSize,
    const int tag,
    const label comm,
    const bool parallel
)
:
    globalIndex()
{
    reset(localSize, tag, comm, parallel);
}


inline Foam::globalIndex::globalIndex(const labelUList& offsets)
:
    offsets_(offsets)
{}


inline Foam::globalIndex::globalIndex(labelList&& offsets)
:
    offsets_(std::move(offsets))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::globalIndex::empty() const
{
    return offsets_.empty() || offsets_.last() == 0;
}


inline const Foam::labelList& Foam::globalIndex::offsets() const noexcept
{
    return offsets_;
}


inline Foam::labelList& Foam::globalIndex::offsets() noexcept
{
    return offsets_;
}


inline Foam::label Foam::globalIndex::size() const
{
    return offsets_.empty() ? 0 : offsets_.last();
}


inline void Foam::globalIndex::reset(const label localSize)
{
    reset(localSize, Pstream::msgType(), UPstream::worldComm, true);
}


inline Foam::label Foam::globalIndex::offset(const label proci) const
{
    return offsets_[proci];
}


inline Foam::label Foam::globalIndex::localStart(const label proci) const
{
    return offsets_[proci];
}


inline Foam::label Foam::globalIndex::localStart() const
{
    return localStart(Pstream::myProcNo());
}


inline Foam::label Foam::globalIndex::localSize(const label proci) const
{
    return offsets_[proci+1] - offsets_[proci];
}


inline Foam::label Foam::globalIndex::localSize() const
{
    return localSize(Pstream::myProcNo());
}


inline Foam::label Foam::globalIndex::maxSize() const
{
    // Use out-of-range proci to avoid excluding any processor
    return maxNonLocalSize(-1);
}


inline Foam::label Foam::globalIndex::maxNonLocalSize() const
{
    return maxNonLocalSize(Pstream::myProcNo());
}


inline Foam::labelRange Foam::globalIndex::range(const label proci) const
{
    return labelRange(offsets_[proci], offsets_[proci+1] - offsets_[proci]);
}


inline Foam::labelRange Foam::globalIndex::range() const
{
    return range(Pstream::myProcNo());
}


inline bool Foam::globalIndex::isLocal(const label proci, const label i) const
{
    return i >= offsets_[proci] && i < offsets_[proci+1];
}


inline bool Foam::globalIndex::isLocal(const label i) const
{
    return isLocal(Pstream::myProcNo(), i);
}


inline Foam::label Foam::globalIndex::toGlobal
(
    const label proci,
    const label i
) const
{
    return i + offsets_[proci];
}


inline Foam::label Foam::globalIndex::toGlobal(const label i) const
{
    return toGlobal(Pstream::myProcNo(), i);
}


inline Foam::labelList Foam::globalIndex::toGlobal
(
    const label proci,
    const labelUList& labels
) const
{
    labelList result(labels);
    inplaceToGlobal(proci, result);

    return result;
}


inline Foam::labelList Foam::globalIndex::toGlobal
(
    const labelUList& labels
) const
{
    return toGlobal(Pstream::myProcNo(), labels);
}


inline void Foam::globalIndex::inplaceToGlobal
(
    const label proci,
    labelList& labels
) const
{
    const label off = offsets_[proci];

    for (label& val : labels)
    {
        val += off;
    }
}


inline void Foam::globalIndex::inplaceToGlobal(labelList& labels) const
{
    inplaceToGlobal(Pstream::myProcNo(), labels);
}


inline Foam::label
Foam::globalIndex::toLocal(const label proci, const label i) const
{
    const label locali = i - offsets_[proci];

    if (locali < 0 || i >= offsets_[proci+1])
    {
        FatalErrorInFunction
            << "Global " << i << " does not belong on processor "
            << proci << nl << "Offsets:" << offsets_
            << abort(FatalError);
    }
    return locali;
}


inline Foam::label Foam::globalIndex::toLocal(const label i) const
{
    return toLocal(Pstream::myProcNo(), i);
}


inline Foam::label Foam::globalIndex::whichProcID(const label i) const
{
    if (i < 0 || i >= size())
    {
        FatalErrorInFunction
            << "Global " << i << " does not belong on any processor."
            << " Offsets:" << offsets_
            << abort(FatalError);
    }

    return findLower(offsets_, i+1);
}


// ************************************************************************* //
