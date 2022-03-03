/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "globalIndex.H"
#include "labelRange.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::globalIndex::reportOverflowAndExit
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


Foam::labelList
Foam::globalIndex::calcOffsets
(
    const labelUList& localLens,
    const bool checkOverflow
)
{
    labelList values;

    const label len = localLens.size();

    if (len)
    {
        values.resize(len+1);

        label start = 0;
        for (label i = 0; i < len; ++i)
        {
            values[i] = start;
            start += localLens[i];

            if (checkOverflow && start < values[i])
            {
                reportOverflowAndExit(i, localLens);
            }
        }
        values[len] = start;
    }

    return values;
}


Foam::List<Foam::labelRange>
Foam::globalIndex::calcRanges
(
    const labelUList& localLens,
    const bool checkOverflow
)
{
    List<labelRange> values;

    const label len = localLens.size();

    if (len)
    {
        values.resize(len);

        label start = 0;
        for (label i = 0; i < len; ++i)
        {
            values[i].reset(start, localLens[i]);
            start += localLens[i];

            if (checkOverflow && start < values[i].start())
            {
                reportOverflowAndExit(i, localLens);
            }
        }
    }

    return values;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::globalIndex::globalIndex(Istream& is)
{
    is >> offsets_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::CompactListList<Foam::label>
Foam::globalIndex::bin
(
    const labelUList& offsets,
    const labelUList& globalIds,
    labelList& order,
    DynamicList<label>& validBins
)
{
    sortedOrder(globalIds, order);
    validBins.clear();

    CompactListList<label> bins;

    if (globalIds.size())
    {
        labelList& binOffsets = bins.offsets();
        binOffsets.resize(offsets.size(), Zero);

        labelList& binValues = bins.values();
        binValues = UIndirectList<label>(globalIds, order);

        const label id = binValues[0];
        label proci = findLower(offsets, id+1);

        validBins.append(proci);
        label binSize = 1;

        for (label i = 1; i < order.size(); i++)
        {
            const label id = binValues[i];

            if (id < offsets[proci+1])
            {
                binSize++;
            }
            else
            {
                // Not local. Reset proci
                label oldProci = proci;
                proci = findLower(offsets, id+1);

                // Set offsets
                for (label j = oldProci+1; j < proci; ++j)
                {
                    binOffsets[j] = binOffsets[oldProci]+binSize;
                }
                binOffsets[proci] = i;
                validBins.append(proci);
                binSize = 1;
            }
        }

        for (label j = proci+1; j < binOffsets.size(); ++j)
        {
            binOffsets[j] = binOffsets[proci]+binSize;
        }
    }

    return bins;
}


void Foam::globalIndex::reset
(
    const label localSize,
    const globalIndex::gatherOnly,
    const label comm
)
{
    // Gather sizes (one-sided)
    reset(UPstream::listGatherValues(localSize, comm));
}


void Foam::globalIndex::reset
(
    const label localSize,
    const label comm,
    const bool parallel
)
{
    labelList localLens;

    const label len = Pstream::nProcs(comm);

    if (len)
    {
        if (parallel && UPstream::parRun())
        {
            localLens = UPstream::listGatherValues(localSize, comm);
            Pstream::broadcast(localLens, comm);
        }
        else
        {
            // Non-parallel branch: use localSize on-proc, zero elsewhere

            localLens.resize(len, Zero);
            localLens[Pstream::myProcNo(comm)] = localSize;
        }

        reset(localLens, true);  // checkOverflow = true
    }
    else
    {
        // Nothing to do
        offsets_.clear();
    }
}


void Foam::globalIndex::reset
(
    const labelUList& localLens,
    const bool checkOverflow
)
{
    const label len = localLens.size();

    if (len)
    {
        offsets_.resize_nocopy(len+1);

        label start = 0;
        for (label i = 0; i < len; ++i)
        {
            offsets_[i] = start;
            start += localLens[i];

            if (checkOverflow && start < offsets_[i])
            {
                reportOverflowAndExit(i, localLens);
            }
        }
        offsets_[len] = start;
    }
    else
    {
        offsets_.clear();
    }
}


void Foam::globalIndex::setLocalSize(const label proci, const label len)
{
    if (proci >= 0 && proci+1 < offsets_.size() && len >= 0)
    {
        const label delta = (len - (offsets_[proci+1] - offsets_[proci]));

        // TBD: additional overflow check
        if (delta)
        {
            for (label i = proci+1; i < offsets_.size(); ++i)
            {
                offsets_[i] += delta;
            }
        }
    }
}


Foam::labelList Foam::globalIndex::localSizes() const
{
    labelList values;

    const label len = (offsets_.size() - 1);

    if (len < 1)
    {
        return values;
    }

    values.resize(len);

    for (label proci=0; proci < len; ++proci)
    {
        values[proci] = offsets_[proci+1] - offsets_[proci];
    }

    return values;
}


Foam::List<Foam::labelRange>
Foam::globalIndex::ranges() const
{
    List<labelRange> values;

    const label len = (offsets_.size() - 1);

    if (len < 1)
    {
        return values;
    }

    values.resize(len);

    for (label proci=0; proci < len; ++proci)
    {
        values[proci].reset
        (
            offsets_[proci],
            (offsets_[proci+1] - offsets_[proci])
        );
    }

    return values;
}


Foam::label Foam::globalIndex::maxNonLocalSize(const label proci) const
{
    const label len = (offsets_.size() - 1);

    if (len < 1)
    {
        return 0;
    }

    label maxLen = 0;

    for (label i=0; i < len; ++i)
    {
        if (i != proci)
        {
            const label localLen = (offsets_[i+1] - offsets_[i]);
            maxLen = max(maxLen, localLen);
        }
    }

    return maxLen;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, globalIndex& gi)
{
    return is >> gi.offsets_;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const globalIndex& gi)
{
    return os << gi.offsets_;
}


// ************************************************************************* //
