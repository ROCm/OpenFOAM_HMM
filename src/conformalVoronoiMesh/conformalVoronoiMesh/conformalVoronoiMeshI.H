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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::conformalVoronoiMesh::timeCheck() const
{
    Info<< nl << "--- [ " << runTime_.elapsedCpuTime() << "s, delta "
        << runTime_.cpuTimeIncrement()<< "s ] --- " << endl;
}


inline Foam::scalar Foam::conformalVoronoiMesh::targetCellSize
(
    const point& pt
) const
{
    return cvMeshControls_.defaultCellSize();
}


inline Foam::scalar Foam::conformalVoronoiMesh::pointPairDistance
(
    const point& pt
) const
{
    return targetCellSize(pt)*cvMeshControls_.pointPairDistanceCoeff();
}


inline Foam::scalar Foam::conformalVoronoiMesh::surfaceSearchDistanceSqr
(
    const point& pt
) const
{
    return sqr(targetCellSize(pt)*cvMeshControls_.surfaceSearchDistanceCoeff());
}


inline Foam::scalar Foam::conformalVoronoiMesh::minimumEdgeLength
(
    const point& pt
) const
{
    return targetCellSize(pt)*cvMeshControls_.minimumEdgeLengthCoeff();
}


inline Foam::tensor Foam::conformalVoronoiMesh::requiredAlignment
(
    const point& pt
) const
{
    Info<< "STILL TO IMPLEMENT requiredAlignment." << endl;
    return I;
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


inline const Foam::Time& Foam::conformalVoronoiMesh::time() const
{
    return runTime_;
}


inline const Foam::conformationSurfaces&
Foam::conformalVoronoiMesh::geometryToConformTo() const
{
    return geometryToConformTo_;
}


inline const Foam::cvControls&
Foam::conformalVoronoiMesh::cvMeshControls() const
{
    return cvMeshControls_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
