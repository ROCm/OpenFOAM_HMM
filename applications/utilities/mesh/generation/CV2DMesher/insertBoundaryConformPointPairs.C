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

#include "CV2D.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::CV2D::insertPointPair
(
    Triangulation::Finite_vertices_iterator& vit,
    const point2D& p,
    const label trii
)
{
    if
    (
        !controls_.mirrorPoints
     || !insertMirrorPoint(toPoint2D(vit->point()), p)
    )
    {
        insertPointPair
        (
            tols_.ppDist,
            p,
            toPoint2D(qSurf_.faceNormals()[trii])
        );
    }

    vit = Triangulation::Finite_vertices_iterator
    (
        CGAL::Filter_iterator
        <
            Triangulation::All_vertices_iterator,
            Triangulation::Infinite_tester
        >(finite_vertices_end(), vit.predicate(), vit.base())
    );
}


bool Foam::CV2D::insertPointPairAtIntersection
(
    Triangulation::Finite_vertices_iterator& vit,
    const point2D& defVert,
    const point2D vertices[],
    const scalar maxProtSize2
)
{
    bool found = false;
    point2D interPoint;
    label interTri = -1;
    scalar interDist2 = 0;

    Face_circulator fcStart = incident_faces(vit);
    Face_circulator fc = fcStart;
    label vi = 0;

    do
    {
        if (!is_infinite(fc))
        {
            pointIndexHit pHit = qSurf_.tree().findLine
            (
                toPoint3D(defVert),
                toPoint3D(vertices[vi])
            );

            if (pHit.hit())
            {
                scalar dist2 =
                    magSqr(toPoint2D(pHit.hitPoint()) - vertices[vi]);

                // Check the point is further away than the furthest so far
                if (dist2 > interDist2)
                {
                    scalar mps2 = maxProtSize2;

                    // If this is a boundary triangle reset the tolerance
                    // to avoid finding a hit point very close to a boundary
                    // vertex
                    if (boundaryTriangle(fc))
                    {
                        mps2 = tols_.maxNotchLen2;
                    }

                    if (dist2 > mps2)
                    {
                        found = true;
                        interPoint = toPoint2D(pHit.hitPoint());
                        interTri = pHit.index();
                        interDist2 = dist2;
                    }
                }
            }

            vi++;
        }
    } while (++fc != fcStart);

    if (found)
    {
        insertPointPair(vit, interPoint, interTri);
        return true;
    }
    else
    {
        return false;
    }
}


Foam::label Foam::CV2D::insertBoundaryConformPointPairs
(
    const fileName& fName
)
{
    label nIntersections = 0;

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        vit++
    )
    {
        // Consider only those points part of point-pairs or near boundary
        if (!vit->nearOrOnBoundary())
        {
            continue;
        }

        // Counter-clockwise circulator
        Face_circulator fcStart = incident_faces(vit);
        Face_circulator fc = fcStart;

        bool infinite = false;
        bool changed = false;

        do
        {
            if (is_infinite(fc))
            {
                infinite = true;
                break;
            }
            else if (fc->faceIndex() < Fb::UNCHANGED)
            {
                changed = true;
                break;
            }
        } while (++fc != fcStart);

        // If the dual-cell is connected to the infinite point or none of the
        // faces whose circumcentres it uses have changed ignore
        if (infinite || !changed) continue;

        fc = fcStart;
        label nVerts = 0;

        do
        {
            vertices[nVerts++] = toPoint2D(circumcenter(fc));

            if (nVerts == maxNvert)
            {
                break;
            }
        } while (++fc != fcStart);

        // Check if dual-cell has a large number of faces in which case
        // assumed to be in the far-field and reject
        if (nVerts == maxNvert) continue;

        // Set n+1 vertex to the first vertex for easy circulating
        vertices[nVerts] = vertices[0];

        // Convert triangle vertex to OpenFOAM point
        point2DFromPoint defVert = toPoint2D(vit->point());

        scalar maxProtSize2 = tols_.maxNotchLen2;

        if (vit->internalOrBoundaryPoint())
        {
            // Calculate metrics of the dual-cell
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            // The perimeter of the dual-cell
            scalar perimeter = 0;

            // Twice the area of the dual-cell
            scalar areaT2 = 0;

            for (int vi=0; vi<nVerts; vi++)
            {
                vector2D edge(vertices[vi+1] - vertices[vi]);
                perimeter += mag(edge);
                vector2D otherEdge = defVert - vertices[vi];
                areaT2 += mag(edge.x()*otherEdge.y() - edge.y()*otherEdge.x());
            }

            // If the dual-cell is very small reject refinement
            if (areaT2 < tols_.minEdgeLen2) continue;

            // Estimate the cell width
            scalar cellWidth = areaT2/perimeter;


            // Check dimensions of dual-cell
            /*
            // Quick rejection of dual-cell refinement based on it's perimeter
            if (perimeter < 2*tols_.minCellSize) continue;

            // Also check the area of the cell and reject refinement
            // if it is less than that allowed
            if (areaT2 < tols_.minCellSize2) continue;

            // Estimate the cell width and reject refinement if it is less than
            // that allowed
            if (cellWidth < 0.5*tols_.minEdgeLen) continue;
            */

            if
            (
                perimeter > 2*controls_.minCellSize
             && areaT2 > controls_.minCellSize2
             && cellWidth > 0.5*tols_.minEdgeLen
            )
            {
                maxProtSize2 = 0.25*tols_.maxNotchLen2;
            }
        }

        if (insertPointPairAtIntersection(vit, defVert, vertices, maxProtSize2))
        {
            nIntersections++;
        }
    }

    return nIntersections;
}


void Foam::CV2D::markNearBoundaryPoints()
{
    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        vit++
    )
    {
        if (vit->internalPoint())
        {
            point vert(toPoint3D(vit->point()));

            pointIndexHit pHit =
                qSurf_.tree().findNearest(vert, 4*controls_.minCellSize2);

            if (pHit.hit())
            {
                vit->setNearBoundary();
            }
        }
    }
}


// ************************************************************************* //
