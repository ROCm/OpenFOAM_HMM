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

#include "CV3D.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::CV3D::dualCellSurfaceIntersection
(
    const Triangulation::Finite_vertices_iterator& vit
) const
{
  // I think that this needs to be done with facets, but think it through.
    return false;
}


void Foam::CV3D::insertSurfaceNearestPointPairs()
{
    Info<< "insertSurfaceNearestPointPairs: " << nl << endl;

    label nSurfacePointsEst = number_of_vertices();

    DynamicList<point> nearSurfacePoints(nSurfacePointsEst);
    DynamicList<point> surfacePoints(nSurfacePointsEst);
    DynamicList<label> surfaceTris(nSurfacePointsEst);

    // Local references to surface mesh addressing
    const pointField& localPoints = qSurf_.localPoints();
    const labelListList& edgeFaces = qSurf_.edgeFaces();
    const vectorField& faceNormals = qSurf_.faceNormals();
    const labelListList& faceEdges = qSurf_.faceEdges();

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        vit++
    )
    {
        if (vit->internalPoint())
        {
	    point vert(topoint(vit->point()));

            pointIndexHit pHit = qSurf_.tree().findNearest
            (
                vert,
                4*controls_.minCellSize2
            );

            if (pHit.hit())
            {
                vit->setNearBoundary();

                // Reference to the nearest triangle
                const labelledTri& f = qSurf_[pHit.index()];

                // Find where point is on triangle.
                // Note tolerance needed is relative one
                // (used in comparing normalized [0..1] triangle coordinates).
                label nearType, nearLabel;
                triPointRef
                (
                    localPoints[f[0]],
                    localPoints[f[1]],
                    localPoints[f[2]]
                ).classify(pHit.hitPoint(), 1e-6, nearType, nearLabel);

                // If point is on a edge check if it is an internal feature

                bool internalFeatureEdge = false;

                if (nearType == triPointRef::EDGE)
                {
                    label edgeI = faceEdges[pHit.index()][nearLabel];
                    const labelList& eFaces = edgeFaces[edgeI];

                    if
                    (
                        eFaces.size() == 2
                     && (faceNormals[eFaces[0]] & faceNormals[eFaces[1]])
                       < -0.2
                    )
                    {
                        internalFeatureEdge = true;
                    }
                }

               if (!internalFeatureEdge && dualCellSurfaceIntersection(vit))
                {
                    nearSurfacePoints.append(vert);
                    surfacePoints.append(pHit.hitPoint());
                    surfaceTris.append(pHit.index());
                }
	    }
	}
    }

    // insertPointPairs
    // (
    //     nearSurfacePoints,
    //     surfacePoints,
    //     surfaceTris,
    //     "surfaceNearestIntersections.obj"
    // );
}


// ************************************************************************* //
