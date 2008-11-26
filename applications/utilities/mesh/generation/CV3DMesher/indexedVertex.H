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

Class
    indexedVertex

Description
    An indexed form of CGAL::Triangulation_vertex_base_3<K> used to keep
    track of the vertices in the triangulation.

\*---------------------------------------------------------------------------*/

#ifndef indexedVertex_H
#define indexedVertex_H

#include <CGAL/Triangulation_3.h>
#include "List.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace CGAL
{

/*---------------------------------------------------------------------------*\
                       Class indexedVertex Declaration
\*---------------------------------------------------------------------------*/

template<class Gt, class Vb=CGAL::Triangulation_vertex_base_3<Gt> >
class indexedVertex
:
    public Vb
{
    // Private data

        //- The index for this triangle vertex
        int index_;

        //- type of pair-point :
        //  NEAR_BOUNDARY_POINT : internal near boundary point.
        //  INTERNAL_POINT      : internal point.
        //  FAR_POINT           : far-point.
        //  >= 0                : part of point-pair. Index of other point.
        //                        Lowest numbered is inside one (master).
        int type_;

        Foam::List<Foam::vector> alignmentDirections_;


public:

    enum pointTypes
    {
        NEAR_BOUNDARY_POINT = -4,
        INTERNAL_POINT      = -3,
        MIRROR_POINT        = -2,
        FAR_POINT           = -1
    };

    typedef typename Vb::Vertex_handle      Vertex_handle;
    typedef typename Vb::Cell_handle        Cell_handle;
    typedef typename Vb::Point              Point;

    template<typename TDS2>
    struct Rebind_TDS
    {
        typedef typename Vb::template Rebind_TDS<TDS2>::Other   Vb2;
        typedef indexedVertex<Gt,Vb2>                           Other;
    };

    indexedVertex()
    :
        Vb(),
        index_(INTERNAL_POINT),
        type_(INTERNAL_POINT),
        alignmentDirections_(0)
    {}

    indexedVertex(const Point& p)
    :
        Vb(p),
        index_(INTERNAL_POINT),
        type_(INTERNAL_POINT),
        alignmentDirections_(0)
    {}

    indexedVertex(const Point& p, Cell_handle f)
    :
        Vb(f, p),
        index_(INTERNAL_POINT),
        type_(INTERNAL_POINT),
        alignmentDirections_(0)
    {}

    indexedVertex(Cell_handle f)
    :
        Vb(f),
        index_(INTERNAL_POINT),
        type_(INTERNAL_POINT),
        alignmentDirections_(0)
    {}


    int& index()
    {
        return index_;
    }

    int index() const
    {
        return index_;
    }


    int& type()
    {
        return type_;
    }

    int type() const
    {
        return type_;
    }

    Foam::List<Foam::vector>& alignmentDirections()
    {
        return alignmentDirections_;
    }

    const Foam::List<Foam::vector>& alignmentDirections() const
    {
        return alignmentDirections_;
    }

    inline bool uninitialised() const
    {
        return type_ == INTERNAL_POINT && index_ == INTERNAL_POINT;
    }

    //- Is point a far-point
    inline bool farPoint() const
    {
        return type_ == FAR_POINT;
    }

    //- Is point internal, i.e. not on boundary
    inline bool internalPoint() const
    {
        return type_ <= INTERNAL_POINT;
    }

    //- Is point internal and near the boundary
    inline bool nearBoundary() const
    {
        return type_ == NEAR_BOUNDARY_POINT;
    }

    //- Set the point to be near the boundary
    inline void setNearBoundary()
    {
        type_ = NEAR_BOUNDARY_POINT;
    }

    //- Is point a mirror point
    inline bool mirrorPoint() const
    {
        return type_ == MIRROR_POINT;
    }

    //- Either master or slave of pointPair.
    inline bool pairPoint() const
    {
        return type_ >= 0;
    }

    //- Master of a pointPair is the lowest numbered one.
    inline bool ppMaster() const
    {
        if (type_ > index_)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    //- Slave of a pointPair is the highest numbered one.
    inline bool ppSlave() const
    {
        if (type_ >= 0 && type_ < index_)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    //- Either original internal point or master of pointPair.
    inline bool internalOrBoundaryPoint() const
    {
        return internalPoint() || ppMaster();
    }

    //- Is point near the boundary or part of the boundary definition
    inline bool nearOrOnBoundary() const
    {
        return pairPoint() || mirrorPoint() || nearBoundary();
    }

    //- Do the two given vertices consitute a boundary point-pair
    inline friend bool pointPair
    (
        const indexedVertex& v0,
        const indexedVertex& v1
    )
    {
        return v0.index_ == v1.type_ || v1.index_ == v0.type_;
    }

    //- Do the three given vertices consitute a boundary triangle
    inline friend bool boundaryTriangle
    (
        const indexedVertex& v0,
        const indexedVertex& v1,
        const indexedVertex& v2
    )
    {
        return (v0.pairPoint() && pointPair(v1, v2))
            || (v1.pairPoint() && pointPair(v2, v0))
            || (v2.pairPoint() && pointPair(v0, v1));
    }

    //- Do the three given vertices consitute an outside triangle
    inline friend bool outsideTriangle
    (
        const indexedVertex& v0,
        const indexedVertex& v1,
        const indexedVertex& v2
    )
    {
        return (v0.farPoint() || v0.ppSlave())
            || (v1.farPoint() || v1.ppSlave())
            || (v2.farPoint() || v2.ppSlave());
    }


    // inline void operator=(const Triangulation::Finite_vertices_iterator vit)
    // {
    //     Vb::operator=indexedVertex(vit->point());

    //     this->index_ = vit->index();

    //     this->type_ = vit->type();

    //     this->displacementSum_ = vit->displacementSum();
    // }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace CGAL


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
