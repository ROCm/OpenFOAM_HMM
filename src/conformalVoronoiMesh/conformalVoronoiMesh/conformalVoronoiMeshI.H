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
