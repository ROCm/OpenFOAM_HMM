/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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

#include "buildCGALPolyhedron.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::buildCGALPolyhedron::buildCGALPolyhedron
(
    const Foam::triSurface& surf
)
:
    CGAL::Modifier_base<HalfedgeDS>(),
    surf_(surf)
{}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::buildCGALPolyhedron::~buildCGALPolyhedron()
{}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::buildCGALPolyhedron::operator()
(
    HalfedgeDS& hds
)
{
    typedef HalfedgeDS::Traits     Traits;
    typedef Traits::Point_3        Point;

    // Postcondition: `hds' is a valid polyhedral surface.
    CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> B(hds, false);

    B.begin_surface
    (
        surf_.points().size(),      // n points
        surf_.size(),               // n facets
        2*surf_.edges().size()      // n halfedges
    );

    forAll(surf_.points(), pI)
    {
        const Foam::point& p = surf_.points()[pI];

        B.add_vertex(Point(p.x(), p.y(), p.z()));
    }

    forAll(surf_, fI)
    {
        B.begin_facet();

        B.add_vertex_to_facet(surf_[fI][0]);
        B.add_vertex_to_facet(surf_[fI][1]);
        B.add_vertex_to_facet(surf_[fI][2]);

        B.end_facet();
    }

    B.end_surface();
}


// ************************************************************************* //
