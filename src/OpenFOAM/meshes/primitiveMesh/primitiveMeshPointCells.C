/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "primitiveMesh.H"
#include "cell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void primitiveMesh::calcPointCells() const
{
    // Loop through cells and mark up points

    if (debug)
    {
        Pout<< "primitiveMesh::calcPointCells() : "
            << "calculating pointCells"
            << endl;
    }

    // It is an error to attempt to recalculate pointCells
    // if the pointer is already set
    if (pcPtr_)
    {
        FatalErrorIn("primitiveMesh::calcPointCells() const")
            << "pointCells already calculated"
            << abort(FatalError);
    }
    else
    {
        const cellList& cf = cells();

        // Count number of cells per point

        labelList npc(nPoints(), 0);

        forAll (cf, cellI)
        {
            const labelList curPoints = cf[cellI].labels(faces());

            forAll (curPoints, pointI)
            {
                label ptI = curPoints[pointI];

                npc[ptI]++;
            }
        }


        // Size and fill cells per point

        pcPtr_ = new labelListList(npc.size());
        labelListList& pointCellAddr = *pcPtr_;

        forAll (pointCellAddr, pointI)
        {
            pointCellAddr[pointI].setSize(npc[pointI]);
        }
        npc = 0;


        forAll (cf, cellI)
        {
            const labelList curPoints = cf[cellI].labels(faces());

            forAll (curPoints, pointI)
            {
                label ptI = curPoints[pointI];

                pointCellAddr[ptI][npc[ptI]++] = cellI;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const labelListList& primitiveMesh::pointCells() const
{
    if (!pcPtr_)
    {
        calcPointCells();
    }

    return *pcPtr_;
}


const labelList& primitiveMesh::pointCells(const label pointI) const
{
    if (hasPointCells())
    {
        return pointCells()[pointI];
    }
    else
    {
        const labelList& own = faceOwner();
        const labelList& nei = faceNeighbour();
        const labelList& pFaces = pointFaces()[pointI];

        labels_.size() = allocSize_;

        label n = 0;

        forAll(pFaces, i)
        {
            const label faceI = pFaces[i];

            // Append owner
            if (n == allocSize_)
            {
                labels_.size() = n;
                allocSize_ = allocSize_*2 + 1;
                labels_.setSize(allocSize_);
            }
            labels_[n++] = own[faceI];

            // Append neighbour
            if (faceI < nInternalFaces())
            {
                if (n == allocSize_)
                {
                    labels_.size() = n;
                    allocSize_ = allocSize_*2 + 1;
                    labels_.setSize(allocSize_);
                }
                labels_[n++] = nei[faceI];
            }
        }
        labels_.size() = n;


        // Filter duplicates
        sort(labels_);

        n = 1;

        for (label i = 1; i < labels_.size(); i++)
        {
            if (labels_[i] != labels_[i-1])
            {
                labels_[n++] = labels_[i];
            }
        }

        labels_.size() = n;

        return labels_;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
