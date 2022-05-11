/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "mapDistributeBase.H"
#include "bitSet.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::mapDistributeBase::renumberMap
(
    labelListList& mapElements,
    const labelUList& oldToNew,
    const bool hasFlip
)
{
    label maxIndex = -1;

    // Transcribe the map
    if (hasFlip)
    {
        for (labelList& map : mapElements)
        {
            for (label& val : map)
            {
                // Unflip indexed value
                const label index = oldToNew[mag(val)-1];

                if (index >= 0)   // Not certain this check is needed
                {
                    maxIndex = max(maxIndex, index);

                    // Retain flip information from original
                    val = (val < 0 ? (-index-1) : (index+1));
                }
            }
        }
    }
    else
    {
        for (labelList& map : mapElements)
        {
            for (label& val : map)
            {
                // Get indexed value (no flipping)

                const label index = oldToNew[val];

                if (index >= 0)   // Not certain this check is needed
                {
                    maxIndex = max(maxIndex, index);
                    val = index;
                }
            }
        }
    }

    return (maxIndex+1);
}


void Foam::mapDistributeBase::renumberVisitOrder
(
    const labelUList& origElements,
    labelList& oldToNew,
    labelListList& maps,
    const bool hasFlip
)
{
    // Both oldToNew and maps refer to compacted numbers in simple
    // ascending order, but we want to recover the original walk order.

    // CAUTION:
    // The following is ill-defined (ie, really bad idea) if the original
    // elements contained duplicates!

    // Inverse mapping:
    //   Original id -> compact id -> walked id

    labelList compactToWalkOrder(origElements.size(), -1);

    forAll(origElements, walkIndex)
    {
        const label origIndex = origElements[walkIndex];
        const label compactIndex = oldToNew[origIndex];

        if (compactIndex >= origElements.size())
        {
            FatalErrorInFunction
                << "Compact index: " << compactIndex
                << " is not represented in the original ("
                << origElements.size()
                << ") elements - indicates an addressing problem" << nl
                << Foam::abort(FatalError);
        }
        else if (compactIndex >= 0)
        {
            compactToWalkOrder[compactIndex] = walkIndex;
            oldToNew[origIndex] = walkIndex;
        }
    }

    renumberMap(maps, compactToWalkOrder, hasFlip);
}


// ************************************************************************* //
