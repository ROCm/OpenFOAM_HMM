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

    As a special exception, you have permission to link this program with the
    CGAL library and distribute executables, as long as you follow the
    requirements of the GNU GPL in regard to all of the software in the
    executable aside from CGAL.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::scalar Foam::conformalVoronoiMesh::targetCellSize
(
    const Foam::point& pt
) const
{
    return cellSizeControl().cellSize(pt);
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
        // There are no internalOrBoundaryPoints available, determine
        // size from scratch

        // Geometric mean
        return sqrt
        (
            targetCellSize(topoint(vA->point()))
           *targetCellSize(topoint(vB->point()))
        );
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
    return targetCellSize(pt)*cvMeshControls().pointPairDistanceCoeff();
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
    return
        sqr
        (
            targetCellSize(pt)
           *cvMeshControls().featurePointExclusionDistanceCoeff()
        );
}


inline Foam::scalar Foam::conformalVoronoiMesh::featureEdgeExclusionDistanceSqr
(
    const Foam::point& pt
) const
{
    return
        sqr
        (
            targetCellSize(pt)
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
    return insertPoint(toPoint(p), type);
}


inline Foam::label Foam::conformalVoronoiMesh::insertPoint
(
    const Point& P,
    const label type
)
{
    uint nVert = number_of_vertices();

    Vertex_handle vh = insert(P);

    if (nVert == number_of_vertices())
    {
        Pout << "Failed to insert point " << topoint(P) << endl;
    }
    else
    {
        vh->index() = nVert;
        vh->type() = type;
    }

    return vh->index();
}


inline void Foam::conformalVoronoiMesh::insertPoint
(
    const Foam::point& p,
    const label index,
    const label type
)
{
    insertPoint(toPoint(p), index, type);
}


inline void Foam::conformalVoronoiMesh::insertPoint
(
    const Point& P,
    const label index,
    const label type
)
{
    uint nVert = number_of_vertices();

    Vertex_handle vh = insert(P);

    if (nVert == number_of_vertices())
    {
        Pout << "Failed to insert point " << topoint(P) << endl;
    }
    else
    {
        vh->index() = index;
        vh->type() = type;
    }
}


inline void Foam::conformalVoronoiMesh::insertPointPair
(
    const scalar ppDist,
    const Foam::point& surfPt,
    const vector& n
)
{
    Pout<< "insertPointPair is depricated, "
        << " it does not check parallel insertion." << endl;

    vector ppDistn = ppDist*n;

    label master = insertPoint
    (
        surfPt - ppDistn,
        number_of_vertices() + 1
    );

    insertPoint(surfPt + ppDistn, master);
}


inline void Foam::conformalVoronoiMesh::createPointPair
(
    const scalar ppDist,
    const Foam::point& surfPt,
    const vector& n,
    DynamicList<Foam::point>& pts,
    DynamicList<label>& indices,
    DynamicList<label>& types
)
{
    vector ppDistn = ppDist*n;

    // Result when the points are eventually inserted.
    // Add number_of_vertices() at insertion of first vertex to all numbers:
    // pt               index type
    // surfPt - ppDistn 0     1
    // surfPt + ppDistn 1     0

    // Master, index = number_of_vertices(), type = number_of_vertices() + 1
    pts.append(surfPt - ppDistn);
    indices.append(0);
    types.append(1);

    // Slave, index = number_of_vertices()(new) + 1, type = index of master.
    // This will be inserted after the master, so number_of_vertices() will
    // have increased
    pts.append(surfPt + ppDistn);
    indices.append(0);
    types.append(-1);
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


inline Foam::List<Foam::label> Foam::conformalVoronoiMesh::processorsAttached
(
    const Delaunay::Finite_facets_iterator& fit
) const
{
    DynamicList<label> procsAttached(8);

    const Cell_handle c1(fit->first);
    const int oppositeVertex = fit->second;
    const Cell_handle c2(c1->neighbor(oppositeVertex));

    FixedList<label, 4> c1Procs(c1->processorsAttached());

    FixedList<label, 4> c2Procs(c2->processorsAttached());

    forAll(c1Procs, aPI)
    {
        if (findIndex(procsAttached, c1Procs[aPI] == -1))
        {
            procsAttached.append(c1Procs[aPI]);
        }

        if (findIndex(procsAttached, c2Procs[aPI] == -1))
        {
            procsAttached.append(c2Procs[aPI]);
        }
    }

    return List<label>(procsAttached);
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


inline const Foam::searchableSurfaces&
Foam::conformalVoronoiMesh::allGeometry() const
{
    return allGeometry_;
}


inline const Foam::conformationSurfaces&
Foam::conformalVoronoiMesh::geometryToConformTo() const
{
    return geometryToConformTo_;
}


inline const Foam::backgroundMeshDecomposition&
Foam::conformalVoronoiMesh::decomposition() const
{
    if (!Pstream::parRun())
    {
        FatalErrorIn
        (
            "inline const Foam::backgroundMeshDecomposition& "
            "Foam::conformalVoronoiMesh::decomposition() const"
        )
            << "The backgroundMeshDecomposition cannot be asked for in serial."
            << exit(FatalError) << endl;
    }

    return decomposition_();
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
