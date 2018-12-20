/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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
    for (label proci = 0; proci < Pstream::nProcs(comm); ++proci)
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


void Foam::globalIndex::reset(const label localSize)
{
    offsets_.resize(Pstream::nProcs()+1);

    labelList localSizes(Pstream::nProcs(), Zero);
    localSizes[Pstream::myProcNo()] = localSize;

    Pstream::gatherList(localSizes, Pstream::msgType());
    Pstream::scatterList(localSizes, Pstream::msgType());

    label offset = 0;
    offsets_[0] = 0;
    for (label proci = 0; proci < Pstream::nProcs(); ++proci)
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
