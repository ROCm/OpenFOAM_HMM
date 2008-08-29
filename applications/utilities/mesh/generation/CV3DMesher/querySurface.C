/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "querySurface.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::querySurface::querySurface
(
    const fileName& surfaceFileName,
    const scalar& includedAngle
)
:
    triSurface(surfaceFileName),
    rndGen_(12345),
    bb_(localPoints()),
    tree_
    (
        treeDataTriSurface(*this),
        bb_.extend(rndGen_, 1e-3), // slightly randomize bb
        8,                         // maxLevel
        4, //10,                   // leafsize
        10.0 //3.0                 // duplicity
    ),
    sFeat_(*this, includedAngle)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::querySurface::~querySurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelListList Foam::querySurface::featurePointFeatureEdges() const
{
    const labelList& featPoints = features().featurePoints();

    const labelList& featEdges = features().featureEdges();

    const edgeList& edges = this->edges();

    List<DynamicList<label> > tempFeatPointFeatEdges(featPoints.size());

    forAll(featPoints, pI)
    {
        label fP = featPoints[pI];

        forAll(featEdges, eI)
        {
            label fE = featEdges[eI];

            const edge& e(edges[fE]);

            if (e.start() == fP || e.end() == fP)
            {
                tempFeatPointFeatEdges[pI].append(fE);
            }
        }
    }

    labelListList featPointFeatEdges(tempFeatPointFeatEdges.size());

    forAll(featPointFeatEdges, fPFE)
    {
        featPointFeatEdges[fPFE].transfer(tempFeatPointFeatEdges[fPFE].shrink());
    }

    return featPointFeatEdges;
}

Foam::indexedOctree<Foam::treeDataTriSurface>::volumeType
Foam::querySurface::insideOutside
(
    const scalar searchSpan2,
    const point& pt
) const
{
    if (!bb_.contains(pt))
    {
        return indexedOctree<treeDataTriSurface>::OUTSIDE;
    }

    pointIndexHit pHit = tree_.findNearest(pt, searchSpan2);

    if (!pHit.hit())
    {
        return tree_.getVolumeType(pt);
    }
    else
    {
        return indexedOctree<treeDataTriSurface>::MIXED;
    }
}


// Check if point is inside surface
bool Foam::querySurface::inside(const point& pt) const
{
    if (!bb_.contains(pt))
    {
        return false;
    }

    return
    (
        tree_.getVolumeType(pt) == indexedOctree<treeDataTriSurface>::INSIDE
    );
}


// Check if point is outside surface
bool Foam::querySurface::outside(const point& pt) const
{
    if (!bb_.contains(pt))
    {
        return true;
    }

    return
    (
        tree_.getVolumeType(pt) == indexedOctree<treeDataTriSurface>::OUTSIDE
    );
}


// Check if point is inside surface by at least dist2
bool Foam::querySurface::wellInside(const point& pt, const scalar dist2) const
{
    if (!bb_.contains(pt))
    {
        return false;
    }

    pointIndexHit pHit = tree_.findNearest(pt, dist2);

    if (pHit.hit())
    {
        return false;
    }
    else
    {
        return
            tree_.getVolumeType(pt)
         == indexedOctree<treeDataTriSurface>::INSIDE;
    }
}


// Check if point is outside surface by at least dist2
bool Foam::querySurface::wellOutside(const point& pt, const scalar dist2) const
{
    if (!bb_.contains(pt))
    {
        return true;
    }

    pointIndexHit pHit = tree_.findNearest(pt, dist2);

    if (pHit.hit())
    {
        return false;
    }
    else
    {
        return
            tree_.getVolumeType(pt)
         == indexedOctree<treeDataTriSurface>::OUTSIDE;
    }
}


bool Foam::querySurface::featurePoint(const label ptI) const
{
    if (findIndex(features().featurePoints(), ptI) >= 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::querySurface::writeTreeOBJ() const
{
    OFstream str("tree.obj");
    label vertI = 0;

    const List<indexedOctree<treeDataTriSurface>::node>& nodes = tree_.nodes();

    forAll(nodes, nodeI)
    {
        const indexedOctree<treeDataTriSurface>::node& nod = nodes[nodeI];

        const treeBoundBox& bb = nod.bb_;

        const pointField points(bb.points());

        label startVertI = vertI;

        forAll(points, i)
        {
            meshTools::writeOBJ(str, points[i]);
            vertI++;
        }

        const edgeList edges(treeBoundBox::edges);

        forAll(edges, i)
        {
            const edge& e = edges[i];

            str << "l " << e[0]+startVertI+1 << ' ' << e[1]+startVertI+1
                << nl;
        }
    }
}


// ************************************************************************* //
