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

Class
    indexedVertex

Description
    An indexed form of CGAL::Triangulation_vertex_base_3<K> used to keep
    track of the Delaunay vertices in the tessellation.

SourceFiles
    indexedVertexI.H
    indexedVertex.C

\*---------------------------------------------------------------------------*/

#ifndef indexedVertex_H
#define indexedVertex_H

#include <CGAL/Triangulation_3.h>
#include "tensor.H"
#include "InfoProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace CGAL
{
template<class Gt, class Vb> class indexedVertex;
}

namespace Foam
{
class Ostream;
template<class Gt, class Vb> Ostream& operator<<
(
    Ostream&,
    const Foam::InfoProxy<CGAL::indexedVertex<Gt, Vb> >&
);
}


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
        //  vtNearBoundary       : internal near boundary point.
        //  vtInternal           : internal point.
        //  vtFar                : far-point.
        //  > vtFar, < 0         : referred point from processor -(type_ + 1)
        //  >= 0                 : part of point-pair. Index of other point.
        //                         Lowest numbered is inside one (master).
        int type_;

        //- Required alignment of the dual cell of this vertex
        Foam::tensor alignment_;

        //- Target size of the dual cell of this vertex
        Foam::scalar targetCellSize_;

        //- Specify whether the vertex is fixed or movable.
        bool vertexFixed_;


public:

    enum vertexTypes
    {
        vtNearBoundary  = INT_MIN,
        vtInternal      = INT_MIN + 1,
        vtFar           = INT_MIN + 2
    };

    enum vertexMotion
    {
        fixed   = 0,
        movable = 1
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


    // Constructors

        inline indexedVertex();

        inline indexedVertex(const Point& p);

        inline indexedVertex(const Point& p, int index, int type);

        inline indexedVertex(const Point& p, Cell_handle f);

        inline indexedVertex(Cell_handle f);


    // Member Functions

        inline int& index();

        inline int index() const;

        inline int& type();

        inline int type() const;

        inline Foam::tensor& alignment();

        inline const Foam::tensor& alignment() const;

        inline Foam::scalar& targetCellSize();

        inline Foam::scalar targetCellSize() const;

        inline bool uninitialised() const;

        //- Is point a far-point
        inline bool farPoint() const;

        //- Is point internal, i.e. not on boundary
        inline bool internalPoint() const;

        //- Is point internal, i.e. not on boundary, external query.
        inline static bool internalPoint(int type);

        // is this a referred vertex
        inline bool referred() const;

        // is this a referred internal or boundary vertex
        inline bool referredInternalOrBoundaryPoint() const;

        // is this a referred external (pair slave) vertex
        inline bool referredExternal() const;

        // Is this a "real" point on this processor, i.e. is it internal or part
        // of the boundary description, and not a "far" or "referred" point
        inline bool real() const;

        // For referred vertices, what is the original processor index
        inline int procIndex() const;

        inline static int encodeProcIndex(int procI);

        //- Set the point to be internal
        inline void setInternal();

        //- Is point internal and near the boundary
        inline bool nearBoundary() const;

        //- Set the point to be near the boundary
        inline void setNearBoundary();

        //- Either master or slave of pointPair.
        inline bool pairPoint() const;

        //- Master of a pointPair is the lowest numbered one.
        inline bool ppMaster() const;

        //- Master of a pointPair is the lowest numbered one, external query.
        inline static bool ppMaster(int index, int type);

        //- Slave of a pointPair is the highest numbered one.
        inline bool ppSlave() const;

        //- Either original internal point or master of pointPair.
        inline bool internalOrBoundaryPoint() const;

        //- Either original internal point or master of pointPair.
        //  External query.
        inline static bool internalOrBoundaryPoint(int index, int type);

        //- Is point near the boundary or part of the boundary definition
        inline bool nearOrOnBoundary() const;

        //- Either a real or referred internal or boundary point
        inline bool anyInternalOrBoundaryPoint() const;

        //- Is the vertex fixed or movable
        inline bool isVertexFixed() const;

        //- Fix the vertex so that it can't be moved
        inline void setVertexFixed() const;

    // inline void operator=(const Delaunay::Finite_vertices_iterator vit)
    // {
    //     Vb::operator=indexedVertex(vit->point());

    //     this->index_ = vit->index();

    //     this->type_ = vit->type();
    // }



    // Info

        //- Return info proxy.
        //  Used to print indexedVertex information to a stream
        Foam::InfoProxy<indexedVertex<Gt, Vb> > info() const
        {
            return *this;
        }

        friend Foam::Ostream& Foam::operator<< <Gt, Vb>
        (
            Foam::Ostream&,
            const Foam::InfoProxy<indexedVertex<Gt, Vb> >&
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace CGAL

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "indexedVertexI.H"

#ifdef NoRepository
#   include "indexedVertex.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
