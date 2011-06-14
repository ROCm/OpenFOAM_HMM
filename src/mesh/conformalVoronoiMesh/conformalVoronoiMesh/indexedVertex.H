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
    indexedVertex

Description
    An indexed form of CGAL::Triangulation_vertex_base_3<K> used to keep
    track of the Delaunay vertices in the tessellation.

\*---------------------------------------------------------------------------*/

#ifndef indexedVertex_H
#define indexedVertex_H

#include <CGAL/Triangulation_3.h>
#include "tensor.H"

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

        //- The index for this Delaunay vertex.  For referred vertices, the
        //  index is negative for vertices that are the outer (slave) of point
        //  pairs
        int index_;

        //- type of pair-point   :
        //  ptNearBoundaryPoint  : internal near boundary point.
        //  ptInternalPoint      : internal point.
        //  ptFarPoint           : far-point.
        //  > ptFarPoint, < 0    : referred point from processor -(type_ + 1)
        //  >= 0                 : part of point-pair. Index of other point.
        //                         Lowest numbered is inside one (master).
        int type_;

        // Required alignment of the dual cell of this vertex
        Foam::tensor alignment_;

        // Target size of the dual cell of this vertex
        Foam::scalar targetCellSize_;


public:

    enum pointTypes
    {
        ptNearBoundaryPoint  = INT_MIN,
        ptInternalPoint      = INT_MIN + 1,
        ptFarPoint           = INT_MIN + 2
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
        index_(ptInternalPoint),
        type_(ptInternalPoint),
        alignment_(),
        targetCellSize_(0.0)
    {}


    indexedVertex(const Point& p)
    :
        Vb(p),
        index_(ptInternalPoint),
        type_(ptInternalPoint),
        alignment_(),
        targetCellSize_(0.0)
    {}


    indexedVertex(const Point& p, int index, int type)
    :
        Vb(p),
        index_(index),
        type_(type),
        alignment_(),
        targetCellSize_(0.0)
    {}


    indexedVertex(const Point& p, Cell_handle f)
    :
        Vb(f, p),
        index_(ptInternalPoint),
        type_(ptInternalPoint),
        alignment_(),
        targetCellSize_(0.0)
    {}


    indexedVertex(Cell_handle f)
    :
        Vb(f),
        index_(ptInternalPoint),
        type_(ptInternalPoint),
        alignment_(),
        targetCellSize_(0.0)
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


    inline Foam::tensor& alignment()
    {
        return alignment_;
    }


    inline const Foam::tensor& alignment() const
    {
        return alignment_;
    }


    inline Foam::scalar& targetCellSize()
    {
        return targetCellSize_;
    }


    inline Foam::scalar targetCellSize() const
    {
        return targetCellSize_;
    }


    inline bool uninitialised() const
    {
        return type_ == ptInternalPoint && index_ == ptInternalPoint;
    }


    //- Is point a far-point
    inline bool farPoint() const
    {
        return type_ == ptFarPoint;
    }


    //- Is point internal, i.e. not on boundary
    inline bool internalPoint() const
    {
        return internalPoint(type_);
    }


    //- Is point internal, i.e. not on boundary, external query.
    inline static bool internalPoint(int type)
    {
        return type <= ptInternalPoint;
    }


    // is this a referred vertex
    inline bool referred() const
    {
        return (type_ < 0 && type_ > ptFarPoint);
    }


    // is this a referred internal or boundary vertex
    inline bool referredInternal() const
    {
        return referred() && index_ >= 0 ;
    }


    // is this a referred external (pair slave) vertex
    inline bool referredExternal() const
    {
        return referred() && index_ < 0;
    }


    // is this a "real" point on this processor, i.e. is it internal or part of
    // the boundary description, and not a "far" or "referred" point
    inline bool real() const
    {
        return internalPoint() || pairPoint();
    }


    // For referred vertices, what is the original processor index
    inline int procIndex() const
    {
        if (referred())
        {
            return -(type_ + 1);
        }
        else
        {
            return -1;
        }
    }


    inline static int encodeProcIndex(int procI)
    {
        return -(procI + 1);
    }


    //- Set the point to be internal
    inline void setInternal()
    {
        type_ = ptInternalPoint;
    }


    //- Is point internal and near the boundary
    inline bool nearBoundary() const
    {
        return type_ == ptNearBoundaryPoint;
    }


    //- Set the point to be near the boundary
    inline void setNearBoundary()
    {
        type_ = ptNearBoundaryPoint;
    }


    //- Either master or slave of pointPair.
    inline bool pairPoint() const
    {
        return type_ >= 0;
    }


    //- Master of a pointPair is the lowest numbered one.
    inline bool ppMaster() const
    {
        return ppMaster(index_, type_);
    }


    //- Master of a pointPair is the lowest numbered one, external query.
    inline static bool ppMaster(int index, int type)
    {
        if (index >= 0 && type > index)
        {
            return true;
        }

        return false;
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
        return internalOrBoundaryPoint(index_, type_);
    }


    //- Either original internal point or master of pointPair, external query.
    inline static bool internalOrBoundaryPoint(int index, int type)
    {
        return internalPoint(type) || ppMaster(index, type);
    }

    //- Is point near the boundary or part of the boundary definition
    inline bool nearOrOnBoundary() const
    {
        return pairPoint() || nearBoundary();
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


    // inline void operator=(const Delaunay::Finite_vertices_iterator vit)
    // {
    //     Vb::operator=indexedVertex(vit->point());

    //     this->index_ = vit->index();

    //     this->type_ = vit->type();
    // }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace CGAL


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
