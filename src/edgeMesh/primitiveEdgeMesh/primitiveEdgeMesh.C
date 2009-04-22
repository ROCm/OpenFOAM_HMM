/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "primitiveEdgeMesh.H"
#include "mergePoints.H"
#include "StaticHashTable.H"
#include "demandDrivenData.H"
#include "meshTools.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::primitiveEdgeMesh::calcPointEdges() const
{
    if (pointEdgesPtr_.valid())
    {
        FatalErrorIn("primitiveEdgeMesh::calcPointEdges() const")
            << "pointEdges already calculated." << abort(FatalError);
    }

    pointEdgesPtr_.reset(new labelListList(points_.size()));
    labelListList& pointEdges = pointEdgesPtr_();

    // Count
    labelList nEdgesPerPoint(points_.size(), 0);

    forAll(edges_, edgeI)
    {
        const edge& e = edges_[edgeI];

        nEdgesPerPoint[e[0]]++;
        nEdgesPerPoint[e[1]]++;
    }

    // Size
    forAll(pointEdges, pointI)
    {
        pointEdges[pointI].setSize(nEdgesPerPoint[pointI]);
    }

    // Fill
    nEdgesPerPoint = 0;

    forAll(edges_, edgeI)
    {
        const edge& e = edges_[edgeI];

        label p0 = e[0];
        pointEdges[p0][nEdgesPerPoint[p0]++] = edgeI;
        label p1 = e[1];
        pointEdges[p1][nEdgesPerPoint[p1]++] = edgeI;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// construct from components
Foam::primitiveEdgeMesh::primitiveEdgeMesh
(
    const pointField& points,
    const edgeList& edges
)
:
    points_(points),
    edges_(edges)
{}


// construct as copy
Foam::primitiveEdgeMesh::primitiveEdgeMesh(const primitiveEdgeMesh& em)
:
    points_(em.points_),
    edges_(em.edges_),
    pointEdgesPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::primitiveEdgeMesh::~primitiveEdgeMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::primitiveEdgeMesh::regions(labelList& edgeRegion) const
{
    edgeRegion.setSize(edges_.size());
    edgeRegion = -1;

    label startEdgeI = 0;

    label currentRegion = 0;

    while (true)
    {
        while (startEdgeI < edges_.size() && edgeRegion[startEdgeI] != -1)
        {
            startEdgeI++;
        }

        if (startEdgeI == edges_.size())
        {
            break;
        }

        // Found edge that has not yet been assigned a region.
        // Mark connected region with currentRegion starting at startEdgeI.

        edgeRegion[startEdgeI] = currentRegion;
        labelList edgesToVisit(1, startEdgeI);

        while (edgesToVisit.size())
        {
            // neighbours of current edgesToVisit
            DynamicList<label> newEdgesToVisit(edgesToVisit.size());

            // Mark all point connected edges with current region.
            forAll(edgesToVisit, i)
            {
                label edgeI = edgesToVisit[i];

                // Mark connected edges
                const edge& e = edges_[edgeI];

                forAll(e, fp)
                {
                    const labelList& pEdges = pointEdges()[e[fp]];

                    forAll(pEdges, pEdgeI)
                    {
                        label nbrEdgeI = pEdges[pEdgeI];

                        if (edgeRegion[nbrEdgeI] == -1)
                        {
                            edgeRegion[nbrEdgeI] = currentRegion;
                            newEdgesToVisit.append(nbrEdgeI);
                        }
                    }
                }
            }

            edgesToVisit.transfer(newEdgesToVisit);
        }

        currentRegion++;
    }
    return currentRegion;
}


void Foam::primitiveEdgeMesh::mergePoints(const scalar mergeDist)
{
    pointField newPoints;
    labelList pointMap;

    bool hasMerged = Foam::mergePoints
    (
        points_,
        mergeDist,
        false,
        pointMap,
        newPoints,
        vector::zero
    );

    if (hasMerged)
    {
        pointEdgesPtr_.clear();

        points_.transfer(newPoints);

        // Renumber and make sure e[0] < e[1] (not really nessecary)
        forAll(edges_, edgeI)
        {
            edge& e = edges_[edgeI];

            label p0 = pointMap[e[0]];
            label p1 = pointMap[e[1]];

            if (p0 < p1)
            {
                e[0] = p0;
                e[1] = p1;
            }
            else
            {
                e[0] = p1;
                e[1] = p0;
            }
        }

        // Compact using a hashtable and commutative hash of edge.
        StaticHashTable<label, edge, Hash<edge> > edgeToLabel
        (
            2*edges_.size()
        );

        label newEdgeI = 0;

        forAll(edges_, edgeI)
        {
            const edge& e = edges_[edgeI];

            if (e[0] != e[1])
            {
                if (edgeToLabel.insert(e, newEdgeI))
                {
                    newEdgeI++;
                }
            }
        }

        edges_.setSize(newEdgeI);

        for
        (
            StaticHashTable<label, edge, Hash<edge> >::const_iterator iter =
                edgeToLabel.begin();
            iter != edgeToLabel.end();
            ++iter
        )
        {
            edges_[iter()] = iter.key();
        }
    }
}


void Foam::primitiveEdgeMesh::writeObj
(
    const fileName& fName
) const
{
    Pout<< nl << "Writing points and edges to " << fName << endl;

    OFstream str(fName);

    forAll(points_, p)
    {
        meshTools::writeOBJ(str, points_[p]);
    }

    forAll (edges_, e)
    {
        const edge& ed = edges_[e];

        str<< "l " << ed.start() + 1 << ' ' << ed.end() + 1;

        str<< nl;
    }
}


// ************************************************************************* //
