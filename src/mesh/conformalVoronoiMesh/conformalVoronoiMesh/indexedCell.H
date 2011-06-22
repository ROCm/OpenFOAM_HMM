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

Class
    indexedCell

Description
    An indexed form of CGAL::Triangulation_cell_base_3<K> used to keep
    track of the Delaunay cells (tets) in the tessellation.

\*---------------------------------------------------------------------------*/

#ifndef indexedCell_H
#define indexedCell_H

#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_cell_base_with_circumcenter_3.h>
#include "indexedVertex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace CGAL
{

/*---------------------------------------------------------------------------*\
                          Class indexedCell Declaration
\*---------------------------------------------------------------------------*/

template
<
    class Gt,
    class Cb=CGAL::Triangulation_cell_base_with_circumcenter_3<Gt>
>
class indexedCell
:
    public Cb
{
    // Private data

        //- The index for this Delaunay tetrahedral cell
        int index_;

        //- The number of times that this Delaunay cell has been limited
        //  during filtering
        int filterCount_;


public:

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
        index_(-1),
        filterCount_(0)
    {}


    indexedCell
    (
        Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Vertex_handle v3
    )
    :
        Cb(v0, v1, v2, v3),
        index_(-1),
        filterCount_(0)
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
        index_(-1),
        filterCount_(0)
    {}


    int& cellIndex()
    {
        return index_;
    }


    int cellIndex() const
    {
        return index_;
    }


    inline int& filterCount()
    {
        return filterCount_;
    }


    inline int filterCount() const
    {
        return filterCount_;
    }


    // Is the Delaunay cell real, i.e. any real vertex
    inline int real() const
    {
        return
        (
            this->vertex(0)->real()
         || this->vertex(1)->real()
         || this->vertex(2)->real()
         || this->vertex(3)->real()
        );
    }


    // Is the Delaunay cell part of the final dual mesj, i.e. any vertex form
    // part of the internal or boundary definition
    inline int internalOrBoundaryDualVertex() const
    {
        return
        (
            this->vertex(0)->internalOrBoundaryPoint()
         || this->vertex(1)->internalOrBoundaryPoint()
         || this->vertex(2)->internalOrBoundaryPoint()
         || this->vertex(3)->internalOrBoundaryPoint()
        );
    }


    // Is the Delaunay cell real or referred (or mixed), i.e. all vertices form
    // part of the real or referred internal or boundary definition
    inline int anyInternalOrBoundaryDualVertex() const
    {
        return
        (
            this->vertex(0)->anyInternalOrBoundaryPoint()
         || this->vertex(1)->anyInternalOrBoundaryPoint()
         || this->vertex(2)->anyInternalOrBoundaryPoint()
         || this->vertex(3)->anyInternalOrBoundaryPoint()
        );
    }


    // A dual vertex on the boundary will result from a Delaunay cell with
    // least one Delaunay vertex outside and at least one inside
    inline bool boundaryDualVertex() const
    {
        return
        (
            (
               this->vertex(0)->internalOrBoundaryPoint()
            || this->vertex(1)->internalOrBoundaryPoint()
            || this->vertex(2)->internalOrBoundaryPoint()
            || this->vertex(3)->internalOrBoundaryPoint()
            )
         && (
               !this->vertex(0)->internalOrBoundaryPoint()
            || !this->vertex(1)->internalOrBoundaryPoint()
            || !this->vertex(2)->internalOrBoundaryPoint()
            || !this->vertex(3)->internalOrBoundaryPoint()
            )
        );
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace CGAL

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
