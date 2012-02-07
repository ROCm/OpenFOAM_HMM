/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Typedefs
    CGALTriangulation3Ddefs

Description
    CGAL data structures used for 3D Delaunay meshing.

    Define CGAL_INEXACT to use Exact_predicates_inexact_constructions kernel
    otherwise the more robust but much less efficient
    Exact_predicates_exact_constructions will be used.

\*---------------------------------------------------------------------------*/

#ifndef CGALTriangulation3Ddefs_H
#define CGALTriangulation3Ddefs_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CGAL/Delaunay_triangulation_3.h"

#include "indexedVertex.H"
#include "indexedCell.H"

#ifdef CGAL_INEXACT
    // Fast kernel using a double as the storage type but the triangulation may
    // fail.  Adding robust circumcentre traits
    #include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
    #include <CGAL/Robust_circumcenter_traits_3.h>
    typedef CGAL::Exact_predicates_inexact_constructions_kernel inexactK;
    typedef CGAL::Robust_circumcenter_traits_3<inexactK>  K;
#else
    // Very robust but expensive kernel
    #include "CGAL/Exact_predicates_exact_constructions_kernel.h"
    typedef CGAL::Exact_predicates_exact_constructions_kernel K;
#endif

typedef CGAL::indexedVertex<K>          Vb;
typedef CGAL::indexedCell<K>            Cb;

typedef CGAL::Triangulation_data_structure_3<Vb, Cb>     Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Delaunay;

typedef Delaunay::Vertex_handle    Vertex_handle;
typedef Delaunay::Cell_handle      Cell_handle;
typedef Delaunay::Point            Point;


//- Spatial sort traits to use with a pair of point pointers and an integer.
//  Taken from a post on the CGAL lists: 2010-01/msg00004.html by
//  Sebastien Loriot (Geometry Factory).
template<class Triangulation>
struct Traits_for_spatial_sort
:
    public Triangulation::Geom_traits
{
    typedef typename Triangulation::Geom_traits Gt;

    typedef std::pair<const typename Triangulation::Point*, int> Point_3;

    struct Less_x_3
    {
        bool operator()(const Point_3& p, const Point_3& q) const
        {
            return typename Gt::Less_x_3()(*(p.first), *(q.first));
        }
    };

    struct Less_y_3
    {
        bool operator()(const Point_3& p, const Point_3& q) const
        {
            return typename Gt::Less_y_3()(*(p.first), *(q.first));
        }
    };

    struct Less_z_3
    {
        bool operator()(const Point_3& p, const Point_3& q) const
        {
            return typename Gt::Less_z_3()(*(p.first), *(q.first));
        }
    };

    Less_x_3  less_x_3_object () const {return Less_x_3();}
    Less_y_3  less_y_3_object () const {return Less_y_3();}
    Less_z_3  less_z_3_object () const {return Less_z_3();}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
