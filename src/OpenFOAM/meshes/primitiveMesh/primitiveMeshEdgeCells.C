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
#include "ListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const labelListList& primitiveMesh::edgeCells() const
{
    if (!ecPtr_)
    {
        if (debug)
        {
            Pout<< "primitiveMesh::edgeCells() : calculating edgeCells" << endl;
        }
        // Invert cellEdges
        ecPtr_ = new labelListList(nEdges());
        invertManyToMany(nEdges(), cellEdges(), *ecPtr_);
    }

    return *ecPtr_;
}


const labelList& primitiveMesh::edgeCells(const label edgeI) const
{
    if (hasEdgeCells())
    {
        return edgeCells()[edgeI];
    }
    else
    {
        const labelList& own = faceOwner();
        const labelList& nei = faceNeighbour();

        // edge faces can either return labels_ or reference in edgeLabels.
        labelList labelsCopy;
        if (!hasEdgeFaces())
        {
            labelsCopy = edgeFaces(edgeI);
        }

        const labelList& eFaces =
        (
            hasEdgeFaces()
          ? edgeFaces()[edgeI]
          : labelsCopy
        );

        labels_.size() = allocSize_;

        // labels_ should certainly be big enough for edge cells.
        label n = 0;

        // Do quadratic insertion.
        forAll(eFaces, i)
        {
            label faceI = eFaces[i];

            {
                label ownCellI = own[faceI];

                // Check if not already in labels_
                for (label j = 0; j < n; j++)
                {
                    if (labels_[j] == ownCellI)
                    {
                        ownCellI = -1;
                        break;
                    }
                }

                if (ownCellI != -1)
                {
                    labels_[n++] = ownCellI;
                }
            }

            if (isInternalFace(faceI))
            {
                label neiCellI = nei[faceI];

                for (label j = 0; j < n; j++)
                {
                    if (labels_[j] == neiCellI)
                    {
                        neiCellI = -1;
                        break;
                    }
                }

                if (neiCellI != -1)
                {
                    labels_[n++] = neiCellI;
                }
            }
        }

        labels_.size() = n;

        return labels_;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
