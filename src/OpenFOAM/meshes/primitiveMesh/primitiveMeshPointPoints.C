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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void primitiveMesh::calcPointPoints() const
{
    if (debug)
    {
        Pout<< "primitiveMesh::calcPointPoints() : "
            << "calculating pointPoints"
            << endl;
    }

    // It is an error to attempt to recalculate pointPoints
    // if the pointer is already set
    if (ppPtr_)
    {
        FatalErrorIn("primitiveMesh::calcPointPoints() const")
            << "pointPoints already calculated"
            << abort(FatalError);
    }
    else
    {
        const edgeList& e = edges();
        const labelListList& pe = pointEdges();

        ppPtr_ = new labelListList(pe.size());
        labelListList& pp = *ppPtr_;

        forAll (pe, pointI)
        {
            pp[pointI].setSize(pe[pointI].size());

            forAll (pe[pointI], ppi)
            {
                if (e[pe[pointI][ppi]].start() == pointI)
                {
                    pp[pointI][ppi] = e[pe[pointI][ppi]].end();
                }
                else if (e[pe[pointI][ppi]].end() == pointI)
                {
                    pp[pointI][ppi] = e[pe[pointI][ppi]].start();
                }
                else
                {
                    FatalErrorIn("primitiveMesh::calcPointPoints() const")
                        << "something wrong with edges"
                        << abort(FatalError);
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const labelListList& primitiveMesh::pointPoints() const
{
    if (!ppPtr_)
    {
        calcPointPoints();
    }

    return *ppPtr_;
}


const labelList& primitiveMesh::pointPoints(const label pointI) const
{
    if (hasPointPoints())
    {
        return pointPoints()[pointI];
    }
    else
    {
        const edgeList& edges = this->edges();
        const labelList& pEdges = pointEdges()[pointI];

        labels_.size() = allocSize_;

        if (pEdges.size() > allocSize_)
        {
            // Set size() so memory allocation behaves as normal.
            labels_.clear();
            allocSize_ = pEdges.size();
            labels_.setSize(allocSize_);
        }

        label n = 0;

        forAll(pEdges, i)
        {
            labels_[n++] = edges[pEdges[i]].otherVertex(pointI);
        }

        labels_.size() = n;

        return labels_;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
