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

Class
    indexedCell

Description
    An indexed form of CGAL::Triangulation_cell_base_3<K> used to keep
    track of the Delaunay cells (tets) in the tessellation.

\*---------------------------------------------------------------------------*/

#ifndef indexedCell_H
#define indexedCell_H

#include <CGAL/Triangulation_3.h>
#include "indexedVertex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace CGAL
{

/*---------------------------------------------------------------------------*\
                       Class indexedCell Declaration
\*---------------------------------------------------------------------------*/

template<class Gt, class Cb=CGAL::Triangulation_cell_base_3<Gt> >
class indexedCell
:
    public Cb
{
    // Private data

        //- The index for this tetrahedral cell
        //  -1: tetrahedral and changed and associated data needs updating
        // >=0: index of cells, set by external update algorithm
        int index_;


public:

    enum cellTypes
    {
        ctUnchanged           = 0,
        ctChanged             = -1,
        ctSaveChanged         = -2
    };

    typedef typename Cb::Vertex_handle                   Vertex_handle;
    typedef typename Cb::Cell_handle                     Cell_handle;

    template < typename TDS2 >
    struct Rebind_TDS
    {
        typedef typename Cb::template Rebind_TDS<TDS2>::Other Cb2;
        typedef indexedCell<Gt, Cb2>  Other;
    };


    indexedCell()
    :
        Cb(),
        index_(ctChanged)
    {}

    indexedCell
    (
        Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Vertex_handle v3
    )
    :
        Cb(v0, v1, v2, v3),
        index_(ctChanged)
    {}

    indexedCell
    (
        Vertex_handle v0,
        Vertex_handle v1,
        Vertex_handle v2,
        Vertex_handle v3,
        Cell_handle n0,
        Cell_handle n1,
        Cell_handle n2,
        Cell_handle n3
    )
    :
        Cb(v0, v1, v2, v3, n0, n1, n2, n3),
        index_(ctChanged)
    {}

    void set_vertex(int i, Vertex_handle v)
    {
        index_ = ctChanged;
        Cb::set_vertex(i, v);
    }

    void set_vertices()
    {
        index_ = ctChanged;
        Cb::set_vertices();
    }

    void set_vertices
    (
        Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Vertex_handle v3
    )
    {
        index_ = ctChanged;
        Cb::set_vertices(v0, v1, v2, v3);
    }

    int& cellIndex()
    {
        return index_;
    }

    int cellIndex() const
    {
        return index_;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace CGAL

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
