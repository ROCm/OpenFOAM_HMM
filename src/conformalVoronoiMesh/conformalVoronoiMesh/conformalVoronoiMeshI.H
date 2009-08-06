/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::conformalVoronoiMesh::timeCheck() const
{
    if (cvMeshControls().timeChecks())
    {
        Info<< nl << "--- [ " << runTime_.elapsedCpuTime() << "s, delta "
            << runTime_.cpuTimeIncrement()<< "s ] --- " << endl;
    }
}


inline Foam::scalar Foam::conformalVoronoiMesh::targetCellSize
(
    const point& pt,
    bool isSurfacePoint
) const
{
    return cellSizeControl().cellSize(pt, isSurfacePoint);
}


inline Foam::scalar Foam::conformalVoronoiMesh::pointPairDistance
(
    const point& pt
) const
{
    // Point pair distances are always going to be at the surface, so the
    // targetCellSize can be told to do a quick, surface only check.

    return targetCellSize(pt, true)*cvMeshControls().pointPairDistanceCoeff();
}


inline Foam::scalar Foam::conformalVoronoiMesh::mixedFeaturePointDistance
(
    const point& pt
) const
{
    return
        pointPairDistance(pt)
       *cvMeshControls().mixedFeaturePointPPDistanceCoeff();
}


inline Foam::scalar Foam::conformalVoronoiMesh::featurePointExclusionDistanceSqr
(
    const point& pt
) const
{
    // Exclusion distance tests are always going to be at the surface, so the
    // targetCellSize can be told to do a quick, surface only check.

    return
        sqr
        (
            targetCellSize(pt, true)
           *cvMeshControls().featurePointExclusionDistanceCoeff()
        );
}


inline Foam::scalar Foam::conformalVoronoiMesh::featureEdgeExclusionDistanceSqr
(
    const point& pt
) const
{
    // Exclusion distance tests are always going to be at the surface, so the
    // targetCellSize can be told to do a quick, surface only check.

    return
        sqr
        (
            targetCellSize(pt, true)
           *cvMeshControls().featureEdgeExclusionDistanceCoeff()
        );
}


inline Foam::scalar Foam::conformalVoronoiMesh::surfaceSearchDistanceSqr
(
    const point& pt
) const
{
    return
        sqr(targetCellSize(pt)*cvMeshControls().surfaceSearchDistanceCoeff());
}


inline Foam::scalar Foam::conformalVoronoiMesh::maxSurfaceProtrusion
(
    const point& pt
) const
{
    return targetCellSize(pt)*cvMeshControls().maxSurfaceProtrusionCoeff();
}


inline Foam::scalar Foam::conformalVoronoiMesh::minimumEdgeLength
(
    const point& pt
) const
{
    return targetCellSize(pt)*cvMeshControls().minimumEdgeLengthCoeff();
}


inline Foam::label Foam::conformalVoronoiMesh::insertPoint
(
    const point& p,
    const label type
)
{
    uint nVert = number_of_vertices();

    Vertex_handle vh = insert(toPoint(p));

    if (nVert == number_of_vertices())
    {
        WarningIn("Foam::conformalVoronoiMesh::insertPoint")
            << "Failed to insert point " << p << endl;
    }
    else
    {
        vh->index() = nVert;
        vh->type() = type;
    }

    return vh->index();
}


inline void Foam::conformalVoronoiMesh::insertPoints
(
    const std::vector<Point>& points
)
{
    Info<< nl << "    " << points.size() << " points to insert..." << endl;

    label nVert = number_of_vertices();

    // using the range insert (faster than inserting points one by one)
    insert(points.begin(), points.end());

    Info<< "    " << number_of_vertices() - startOfInternalPoints_
        << " points inserted" << endl;

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->uninitialised())
        {
            vit->index() = nVert++;
        }
    }
}


inline void Foam::conformalVoronoiMesh::insertPointPair
(
    const scalar ppDist,
    const point& surfPt,
    const vector& n
)
{
    vector ppDistn = ppDist*n;

    label master = insertPoint
    (
        surfPt - ppDistn,
        number_of_vertices() + 1
    );

    insertPoint(surfPt + ppDistn, master);
}


inline void Foam::conformalVoronoiMesh::insertVb(const Vb& v, label offset)
{
    const Point& Pt(v.point());

    uint nVert = number_of_vertices();

    Vertex_handle vh = insert(Pt);

    if (nVert == number_of_vertices())
    {
        FatalErrorIn("Foam::conformalVoronoiMesh::insertVb(const Vb& v")
            << "Failed to reinsert Vb at " << topoint(Pt)
            << endl;
    }
    else
    {
        vh->index() = v.index() + offset;

        if (v.pairPoint())
        {
            vh->type() = v.type() + offset;
        }
        else
        {
            vh->type() = v.type();
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#ifdef CGAL_INEXACT

inline Foam::conformalVoronoiMesh::pointFromPoint
Foam::conformalVoronoiMesh::topoint
(
    const Point& P
) const
{
    return reinterpret_cast<pointFromPoint>(P);
}


inline Foam::conformalVoronoiMesh::PointFrompoint
Foam::conformalVoronoiMesh::toPoint
(
    const point& p
) const
{
    return reinterpret_cast<PointFrompoint>(p);
}


#else

inline Foam::conformalVoronoiMesh::pointFromPoint
Foam::conformalVoronoiMesh::topoint
(
    const Point& P
) const
{
    return point
    (
        CGAL::to_double(P.x()),
        CGAL::to_double(P.y()),
        CGAL::to_double(P.z())
    );
}


inline Foam::conformalVoronoiMesh::PointFrompoint
Foam::conformalVoronoiMesh::toPoint
(
    const point& p
) const
{
    return Point(p.x(), p.y(), p.z());
}

#endif


inline Foam::conformalVoronoiMesh::CGALVector
Foam::conformalVoronoiMesh::toCGALVector(const point& pt) const
{
    return CGALVector(pt.x(), pt.y(), pt.z());
}


inline const Foam::Time& Foam::conformalVoronoiMesh::time() const
{
    return runTime_;
}


inline const Foam::conformationSurfaces&
Foam::conformalVoronoiMesh::geometryToConformTo() const
{
    return geometryToConformTo_;
}


inline const Foam::cellSizeControlSurfaces&
Foam::conformalVoronoiMesh::cellSizeControl() const
{
    return cellSizeControl_;
}


inline const Foam::cvControls&
Foam::conformalVoronoiMesh::cvMeshControls() const
{
    return cvMeshControls_;
}


// ************************************************************************* //
