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

#include "globalIndex.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::globalIndex::globalIndex(Istream& is)
{
    is >> offsets_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::globalIndex::bin
(
    const labelUList& offsets,
    const labelUList& globalIds,
    labelList& order,
    CompactListList<label>& bins,
    DynamicList<label>& validBins
)
{
    sortedOrder(globalIds, order);

    bins.m() = UIndirectList<label>(globalIds, order);

    labelList& binOffsets = bins.offsets();
    binOffsets.setSize(offsets.size());
    binOffsets = 0;

    validBins.clear();

    if (globalIds.size())
    {
        const label id = bins.m()[0];
        label proci = findLower(offsets, id+1);

        validBins.append(proci);
        label binSize = 1;

        for (label i = 1; i < order.size(); i++)
        {
            const label id = bins.m()[i];

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
}


void Foam::globalIndex::reset
(
    const label localSize,
    const int tag,
    const label comm,
    const bool parallel
)
{
    offsets_.resize(Pstream::nProcs(comm)+1);

    labelList localSizes(Pstream::nProcs(comm), Zero);
    localSizes[Pstream::myProcNo(comm)] = localSize;

    if (parallel)
    {
        Pstream::gatherList(localSizes, tag, comm);
        Pstream::scatterList(localSizes, tag, comm);
    }

    label offset = 0;
    offsets_[0] = 0;
    for (const int proci : Pstream::allProcs(comm))
    {
        const label oldOffset = offset;
        offset += localSizes[proci];

        if (offset < oldOffset)
        {
            FatalErrorInFunction
                << "Overflow : sum of sizes " << localSizes
                << " exceeds capability of label (" << labelMax
                << "). Please recompile with larger datatype for label."
                << exit(FatalError);
        }
        offsets_[proci+1] = offset;
    }
}


Foam::labelList Foam::globalIndex::sizes() const
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
