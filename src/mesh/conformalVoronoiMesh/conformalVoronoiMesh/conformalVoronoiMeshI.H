/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
     \\/     M anipulation  |
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::scalar Foam::conformalVoronoiMesh::targetCellSize
(
    const Foam::point& pt,
    bool isSurfacePoint
) const
{
    return cellSizeControl().cellSize(pt, isSurfacePoint);
}


inline Foam::scalar Foam::conformalVoronoiMesh::averageCellSize
(
    const Vertex_handle& vA,
    const Vertex_handle& vB
) const
{
    // Arithmetic mean
    // return 0.5*(vA->targetCellSize() + vB->targetCellSize());

    // Geometric mean
    return sqrt(vA->targetCellSize()*vB->targetCellSize());

    // Harmonic mean
    // return
    //     2.0*(vA->targetCellSize()*vB->targetCellSize())
    //    /(vA->targetCellSize() + vB->targetCellSize());
}


inline Foam::scalar Foam::conformalVoronoiMesh::averageAnyCellSize
(
    const Vertex_handle& vA,
    const Vertex_handle& vB
) const
{
    if
    (
        !vA->internalOrBoundaryPoint()
     && !vB->internalOrBoundaryPoint()
    )
    {
        FatalErrorIn
        (
            "Foam::conformalVoronoiMesh::averageAnyCellSize"
            "("
                "const Vertex_handle& vA,"
                "const Vertex_handle& vB"
            ") const"
        )
            << "Requested averageCellSize for two external points"
            << nl
            << exit(FatalError);
    }
    else if (!vB->internalOrBoundaryPoint())
    {
        return vA->targetCellSize();
    }
    else if (!vA->internalOrBoundaryPoint())
    {
        return vB->targetCellSize();
    }

    return averageCellSize(vA, vB);
}


inline Foam::scalar Foam::conformalVoronoiMesh::averageAnyCellSize
(
    const Delaunay::Finite_facets_iterator& fit
) const
{
    // Geometric mean

    scalar sizeProduct = 1;
    label nProducts = 0;

    const Cell_handle c(fit->first);
    const int oppositeVertex = fit->second;

    for (label i = 0; i < 3; i++)
    {
        Vertex_handle v = c->vertex(vertex_triple_index(oppositeVertex, i));

        if (v->internalOrBoundaryPoint())
        {

            sizeProduct *= v->targetCellSize();

            nProducts++;
        }
    }

    if (nProducts < 1)
    {
        // There are no internalOrBoundaryPoints available, determine
        // size from scratch

        for (label i = 0; i < 3; i++)
        {
            Vertex_handle v = c->vertex(vertex_triple_index(oppositeVertex, i));

            sizeProduct *= targetCellSize(topoint(v->point()));
        }

        nProducts = 3;
    }

    return pow(sizeProduct, (1.0/nProducts));
}


inline Foam::scalar Foam::conformalVoronoiMesh::pointPairDistance
(
    const Foam::point& pt
) const
{
    // Point pair distances are always going to be at the surface, so the
    // targetCellSize can be told to do a quick, surface only check.

    return targetCellSize(pt, true)*cvMeshControls().pointPairDistanceCoeff();
}


inline Foam::scalar Foam::conformalVoronoiMesh::mixedFeaturePointDistance
(
    const Foam::point& pt
) const
{
    return
        pointPairDistance(pt)
       *cvMeshControls().mixedFeaturePointPPDistanceCoeff();
}


inline Foam::scalar Foam::conformalVoronoiMesh::featurePointExclusionDistanceSqr
(
    const Foam::point& pt
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
    const Foam::point& pt
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
    const Foam::point& pt
) const
{
    return
        sqr(targetCellSize(pt)*cvMeshControls().surfaceSearchDistanceCoeff());
}


inline Foam::scalar Foam::conformalVoronoiMesh::maxSurfaceProtrusion
(
    const Foam::point& pt
) const
{
    return targetCellSize(pt)*cvMeshControls().maxSurfaceProtrusionCoeff();
}


inline Foam::label Foam::conformalVoronoiMesh::insertPoint
(
    const Foam::point& p,
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
    Info<< "    " << points.size() << " points to insert..." << endl;

    label nVert = number_of_vertices();

    // using the range insert (faster than inserting points one by one)
    insert(points.begin(), points.end());

    Info<< "    " << number_of_vertices() - startOfInternalPoints_
        << " points inserted" << endl;

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
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
    const Foam::point& surfPt,
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
            << nl
            << exit(FatalError);
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


inline bool Foam::conformalVoronoiMesh::isBoundaryDualFace
(
    const Delaunay::Finite_edges_iterator& eit
) const
{
    Cell_handle c = eit->first;
    Vertex_handle vA = c->vertex(eit->second);
    Vertex_handle vB = c->vertex(eit->third);

    // A dual face on the boundary will result from one Dv inside and
    // one outside
    return
    (
        (
            vA->internalOrBoundaryPoint()
         || vB->internalOrBoundaryPoint()
        )
     && (
            !vA->internalOrBoundaryPoint()
         || !vB->internalOrBoundaryPoint()
        )
    );
}


inline bool Foam::conformalVoronoiMesh::isBoundaryDualVertex
(
    const Delaunay::Finite_cells_iterator& cit
) const
{
    // A dual vertex on the boundary will result from a Delaunay cell
    // with least one Delaunay vertex outside and at least one inside

    return
    (
        (
            cit->vertex(0)->internalOrBoundaryPoint()
         || cit->vertex(1)->internalOrBoundaryPoint()
         || cit->vertex(2)->internalOrBoundaryPoint()
         || cit->vertex(3)->internalOrBoundaryPoint()
        )
     && (
            !cit->vertex(0)->internalOrBoundaryPoint()
         || !cit->vertex(1)->internalOrBoundaryPoint()
         || !cit->vertex(2)->internalOrBoundaryPoint()
         || !cit->vertex(3)->internalOrBoundaryPoint()
        )
    );
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
    const Foam::point& p
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
    return Foam::point
    (
        CGAL::to_double(P.x()),
        CGAL::to_double(P.y()),
        CGAL::to_double(P.z())
    );
}


inline Foam::conformalVoronoiMesh::PointFrompoint
Foam::conformalVoronoiMesh::toPoint
(
    const Foam::point& p
) const
{
    return Point(p.x(), p.y(), p.z());
}

#endif


inline Foam::conformalVoronoiMesh::CGALVector
Foam::conformalVoronoiMesh::toCGALVector(const Foam::vector& v) const
{
    return CGALVector(v.x(), v.y(), v.z());
}


inline const Foam::Time& Foam::conformalVoronoiMesh::time() const
{
    return runTime_;
}


inline Foam::Random& Foam::conformalVoronoiMesh::rndGen() const
{
    return rndGen_;
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
