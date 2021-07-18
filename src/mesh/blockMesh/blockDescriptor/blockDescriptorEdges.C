/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "blockDescriptor.H"
#include "lineEdge.H"
#include "lineDivide.H"

// * * * * * * * * * * * * * * Local Data Members  * * * * * * * * * * * * * //

// Warning.
// Ordering of edges needs to be the same as hex cell shape model
static const int hexEdge0[12] = { 0, 3, 7, 4,  0, 1, 5, 4,  0, 1, 2, 3 };
static const int hexEdge1[12] = { 1, 2, 6, 5,  3, 2, 6, 7,  4, 5, 6, 7 };


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

int Foam::blockDescriptor::calcEdgePointsWeights
(
    pointField& edgePoints,
    scalarList& edgeWeights,
    const label start,
    const label end,
    const label nDiv,
    const gradingDescriptors& expand
) const
{
    // Set reference to the list of labels defining the block
    const labelList& blockLabels = blockShape_;

    // Set the edge points/weights
    // The edge is a straight-line if it is not in the list of blockEdges

    for (const blockEdge& cedge : blockEdges_)
    {
        const int cmp = cedge.compare(blockLabels[start], blockLabels[end]);

        if (cmp > 0)
        {
            // Curve has the same orientation

            // Divide the line
            const lineDivide divEdge(cedge, nDiv, expand);

            edgePoints  = divEdge.points();
            edgeWeights = divEdge.lambdaDivisions();

            return 1;  // Found curved-edge: done
        }
        else if (cmp < 0)
        {
            // Curve has the opposite orientation

            // Divide the line
            const lineDivide divEdge(cedge, nDiv, expand.inv());

            const pointField& p = divEdge.points();
            const scalarList& d = divEdge.lambdaDivisions();

            edgePoints.resize(p.size());
            edgeWeights.resize(d.size());

            // Copy in reverse order
            const label pn = (p.size() - 1);
            forAll(p, pi)
            {
                edgePoints[pi]  = p[pn - pi];
                edgeWeights[pi] = 1 - d[pn - pi];
            }

            return 1;  // Found curved-edge: done
        }
    }

    // Not curved-edge: divide the edge as a straight line

    // Get list of points for this block
    const pointField blockPoints(blockShape_.points(vertices_));

    lineDivide divEdge
    (
        blockEdges::lineEdge(blockPoints, start, end),
        nDiv,
        expand
    );

    edgePoints = divEdge.points();
    edgeWeights = divEdge.lambdaDivisions();

    return 0;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

int Foam::blockDescriptor::edgesPointsWeights
(
    pointField (&edgesPoints)[12],
    scalarList (&edgesWeights)[12]
) const
{
    int nCurved = 0;

    for (label edgei = 0; edgei < 12; ++edgei)
    {
        nCurved += calcEdgePointsWeights
        (
            edgesPoints[edgei],
            edgesWeights[edgei],
            hexEdge0[edgei],
            hexEdge1[edgei],

            sizes()[edgei/4],   // 12 edges -> 3 components (x,y,z)
            expand_[edgei]
        );
    }

    return nCurved;
}


bool Foam::blockDescriptor::edgePointsWeights
(
    const label edgei,
    pointField& edgePoints,
    scalarList& edgeWeights,
    const label nDiv,
    const gradingDescriptors& gd
) const
{
    if (edgei < 0 || edgei >= 12)
    {
        FatalErrorInFunction
            << "Edge label " << edgei
            << " out of range 0..11"
            << exit(FatalError);
    }

    const int nCurved = calcEdgePointsWeights
    (
        edgePoints,
        edgeWeights,
        hexEdge0[edgei],
        hexEdge1[edgei],
        nDiv,
        gd
    );

    return nCurved;
}


bool Foam::blockDescriptor::edgePointsWeights
(
    const label edgei,
    pointField& edgePoints,
    scalarList& edgeWeights
) const
{
    if (edgei < 0 || edgei >= 12)
    {
        FatalErrorInFunction
            << "Edge label " << edgei
            << " out of range 0..11"
            << exit(FatalError);
    }

    const int nCurved = calcEdgePointsWeights
    (
        edgePoints,
        edgeWeights,
        hexEdge0[edgei],
        hexEdge1[edgei],
        sizes()[edgei/4],   // 12 edges -> 3 components (x,y,z)
        expand_[edgei]
    );

    return nCurved;
}


// ************************************************************************* //
