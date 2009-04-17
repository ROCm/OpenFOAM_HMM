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

#include "featureEdgeMesh.H"
#include "Random.H"
#include "meshTools.H"
#include "linePointRef.H"
#include "OFstream.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::featureEdgeMesh::featureEdgeMesh(const IOobject& io)
:
    edgeMesh(io),
    normals_(),
    edgeNormals_(edges().size()),
    allEdges_(identity(edges().size()))

{}


Foam::featureEdgeMesh::featureEdgeMesh
(
    const IOobject& io,
    const pointField& pts,
    const edgeList& eds,
    const vectorField& normals
)
:
    edgeMesh(io, pts, eds),
    normals_(normals),
    edgeNormals_(edges().size()),
    allEdges_(identity(edges().size()))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::featureEdgeMesh::~featureEdgeMesh()
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::featureEdgeMesh::nearestFeatureEdge
(
    const pointField& samples,
    scalarField searchDistSqr,
    labelList& edgeLabel,
    pointField& edgePoint,
    labelListList& adjacentNormals
) const
{
    edgeLabel.setSize(samples.size());
    edgePoint.setSize(samples.size());
    adjacentNormals.setSize(samples.size());

    forAll(samples, i)
    {
        const point& sample = samples[i];

        pointIndexHit pHit = edgeTree().findNearest
        (
            sample,
            searchDistSqr[i]
        );

        if (!pHit.hit())
        {
            edgeLabel[i] = -1;
        }
        else
        {
            edgeLabel[i] = allEdges_[pHit.index()];
            edgePoint[i] = pHit.hitPoint();
            adjacentNormals[i] = edgeNormals_[edgeLabel[i]];
        }
    }
}


const Foam::indexedOctree<Foam::treeDataEdge>&
Foam::featureEdgeMesh::edgeTree() const
{
    if (edgeTree_.empty())
    {
        Random rndGen(17301893);

        // Slightly extended bb. Slightly off-centred just so on symmetric
        // geometry there are less face/edge aligned items.
        treeBoundBox bb
        (
            treeBoundBox(points()).extend(rndGen, 1E-4)
        );

        bb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        bb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

        edgeTree_.reset
        (
            new indexedOctree<treeDataEdge>
            (
                treeDataEdge
                (
                    false,          // cachebb
                    edges(),        // edges
                    points(),       // points
                    allEdges_       // selected edges
                ),
                bb,     // bb
                8,      // maxLevel
                10,     // leafsize
                3.0     // duplicity
            )
        );
    }

    return edgeTree_();
}


// ************************************************************************* //
