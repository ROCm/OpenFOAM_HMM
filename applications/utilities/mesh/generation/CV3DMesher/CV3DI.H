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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::label Foam::CV3D::insertPoint
(
    const point& p,
    const label type
)
{
    uint nVert = number_of_vertices();

    Vertex_handle vh = insert(toPoint(p));

    if (nVert == number_of_vertices())
    {
        WarningIn("Foam::CV3D::insertPoint")
            << "Failed to insert point " << p << endl;
    }
    else
    {
        vh->index() = nVert;
        vh->type() = type;
    }

    return vh->index();
}


inline bool Foam::CV3D::insertMirrorPoint
(
    const point& nearSurfPt,
    const point& surfPt
)
{
    point mirrorPoint(2*surfPt - nearSurfPt);

    if (qSurf_.outside(mirrorPoint))
    {
        insertPoint(mirrorPoint, Vb::MIRROR_POINT);
        return true;
    }
    else
    {
        return false;
    }
}


inline void Foam::CV3D::insertPointPair
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


inline void Foam::CV3D::insertVb(const Vb& v)
{
    const Point& Pt(v.point());

    uint nVert = number_of_vertices();

    Vertex_handle vh = insert(Pt);

    if (nVert == number_of_vertices())
    {
        FatalErrorIn("Foam::CV3D::insert(const Vb& v")
        << "Failed to reinsert Vb at " << topoint(Pt)
            << endl;
    }
    else
    {
        vh->index() = v.index();
        vh->type() = v.type();
    }
}


void Foam::CV3D::movePoint(const Vertex_handle& vh, const point& p)
{
    label vi = vh->index();

    // remove(vh);
    // insert(toPoint(p))->index() = vi;

    move_point(vh, toPoint(p));
    vh->index() = vi;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#ifdef CGAL_INEXACT

inline Foam::CV3D::pointFromPoint Foam::CV3D::topoint
(
    const Point& P
) const
{
    return reinterpret_cast<pointFromPoint>(P);
}

inline Foam::CV3D::PointFrompoint Foam::CV3D::toPoint
(
    const point& p
) const
{
    return reinterpret_cast<PointFrompoint>(p);
}

#else

inline Foam::CV3D::pointFromPoint Foam::CV3D::topoint
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

inline Foam::CV3D::PointFrompoint Foam::CV3D::toPoint
(
    const point& p
) const
{
    return Point(p.x(), p.y(), p.z());
}

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
