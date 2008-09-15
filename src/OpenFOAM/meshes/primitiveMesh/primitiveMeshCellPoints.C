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

const labelListList& primitiveMesh::cellPoints() const
{
    if (!cpPtr_)
    {
        if (debug)
        {
            Pout<< "primitiveMesh::cellPoints() : "
                << "calculating cellPoints" << endl;
        }

        // Invert pointCells
        cpPtr_ = new labelListList(nCells());
        invertManyToMany(nCells(), pointCells(), *cpPtr_);
    }
    
    return *cpPtr_;
}


const labelList& primitiveMesh::cellPoints(const label cellI) const
{
    if (hasCellPoints())
    {
        return cellPoints()[cellI];
    }
    else
    {
        const faceList& fcs = faces();
        const labelList& cFaces = cells()[cellI];

        labelSet_.clear();

        forAll(cFaces, i)
        {
            const labelList& f = fcs[cFaces[i]];

            forAll(f, fp)
            {
                labelSet_.insert(f[fp]);    
            }
        }

        labels_.size() = allocSize_;

        if (labelSet_.size() > allocSize_)
        {
            labels_.clear();
            allocSize_ = labelSet_.size();
            labels_.setSize(allocSize_);
        }

        label n = 0;

        forAllConstIter(labelHashSet, labelSet_, iter)
        {
            labels_[n++] = iter.key();
        }

        labels_.size() = n;

        return labels_;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
