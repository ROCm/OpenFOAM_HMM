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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::label Foam::CV2D::insertPoint
(
    const point2D& p,
    const label type
)
{
    uint nVert = number_of_vertices();

    Vertex_handle vh = insert(toPoint(p));

    if (nVert == number_of_vertices())
    {
        WarningIn("Foam::CV2D::insertPoint")
            << "Failed to insert point " << p << endl;
    }
    else
    {
        vh->index() = nVert;
        vh->type() = type;
    }

    return vh->index();
}


inline bool Foam::CV2D::insertMirrorPoint
(
    const point2D& nearSurfPt,
    const point2D& surfPt
)
{
    point2D mirrorPoint(2*surfPt - nearSurfPt);

    if (qSurf_.outside(toPoint3D(mirrorPoint)))
    {
        insertPoint(mirrorPoint, Vb::MIRROR_POINT);
        return true;
    }
    else
    {
        return false;
    }
}


inline void Foam::CV2D::insertPointPair
(
    const scalar ppDist,
    const point2D& surfPt,
    const vector2D& n
)
{
    vector2D ppDistn = ppDist*n;

    label master = insertPoint
    (
        surfPt - ppDistn,
        number_of_vertices() + 1
    );
    
    insertPoint(surfPt + ppDistn, master);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::point2D& Foam::CV2D::toPoint2D(const point& p) const
{
    return reinterpret_cast<const point2D&>(p);
}

inline Foam::point Foam::CV2D::toPoint3D(const point2D& p) const
{
    return point(p.x(), p.y(), z_);
}


#ifdef CGAL_INEXACT

inline Foam::CV2D::point2DFromPoint Foam::CV2D::toPoint2D(const Point& P) const
{
    return reinterpret_cast<point2DFromPoint>(P);
}

inline Foam::CV2D::PointFromPoint2D Foam::CV2D::toPoint(const point2D& p) const
{
    return reinterpret_cast<PointFromPoint2D>(p);
}

#else

inline Foam::CV2D::point2DFromPoint Foam::CV2D::toPoint2D(const Point& P) const
{
    return point2D(CGAL::to_double(P.x()), CGAL::to_double(P.y()));
}

inline Foam::CV2D::PointFromPoint2D Foam::CV2D::toPoint(const point2D& p) const
{
    return Point(p.x(), p.y());
}

#endif


inline Foam::point Foam::CV2D::toPoint3D(const Point& P) const
{
    return point(CGAL::to_double(P.x()), CGAL::to_double(P.y()), z_);
}


inline void Foam::CV2D::movePoint(const Vertex_handle& vh, const point2D& p)
{
    vh->set_point(toPoint(p));
    fast_restore_Delaunay(vh);
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

inline bool Foam::boundaryTriangle(const CV2D::Face_handle fc)
{
    return boundaryTriangle
    (
        *fc->vertex(0),
        *fc->vertex(1),
        *fc->vertex(2)
    );
}

inline bool Foam::outsideTriangle(const CV2D::Face_handle fc)
{
    return outsideTriangle
    (
        *fc->vertex(0),
        *fc->vertex(1),
        *fc->vertex(2)
    );
}


// ************************************************************************* //
