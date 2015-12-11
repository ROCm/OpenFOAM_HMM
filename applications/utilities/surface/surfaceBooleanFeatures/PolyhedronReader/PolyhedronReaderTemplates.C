/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

#include "PolyhedronReader.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class HDS>
Foam::PolyhedronReader::Build_triangle<HDS>::Build_triangle
(
    const triSurface& s
)
:
    s_(s)
{}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class HDS>
void Foam::PolyhedronReader::Build_triangle<HDS>::operator()(HDS& hds)
{
    // Postcondition: hds is a valid polyhedral surface.
    CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);

    B.begin_surface(s_.nPoints(), s_.size());

    typedef typename HDS::Vertex Vertex;
    typedef typename Vertex::Point Point;

    const Foam::pointField& pts = s_.points();
    forAll(pts, i)
    {
        const Foam::point& pt = pts[i];
        B.add_vertex(Point(pt.x(), pt.y(), pt.z()));
    }
    forAll(s_, i)
    {
        const Foam::labelledTri& t = s_[i];
        B.begin_facet();
        B.add_vertex_to_facet(t[0]);
        B.add_vertex_to_facet(t[1]);
        B.add_vertex_to_facet(t[2]);
        B.end_facet();
    }
    B.end_surface();
}


// ************************************************************************* //
