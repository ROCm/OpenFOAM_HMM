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

    As a special exception, you have permission to link this program with the
    CGAL library and distribute executables, as long as you follow the
    requirements of the GNU GPL in regard to all of the software in the
    executable aside from CGAL.

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
#include "List.H"
#include "globalIndex.H"
#include "Pstream.H"
#include "Swap.H"
#include "InfoProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace CGAL
{
template<class Gt, class Cb> class indexedCell;
}

namespace Foam
{
class Ostream;
template<class Gt, class Cb> Ostream& operator<<
(
    Ostream&,
    const Foam::InfoProxy<CGAL::indexedCell<Gt, Cb> >&
);
}


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

        //- The index for this Delaunay tetrahedral cell.  Type information is
        //- also carried:

        //  ctFar : the dual point of this cell does not form part of the
        //          internal or boundary of the dual mesh
        //  >=0   : the (local) index of an internal or boundary dual point,
        //           not on a processor face
        //  < 0 && > ctFar : the (global) index of a dual point on a processor
        //                   face
        int index_;

        //- The number of times that this Delaunay cell has been limited
        //  during filtering
        int filterCount_;


public:

    enum cellTypes
    {
        ctFar = INT_MIN
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
        index_(ctFar),
        filterCount_(0)
    {}


    indexedCell
    (
        Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Vertex_handle v3
    )
    :
        Cb(v0, v1, v2, v3),
        index_(ctFar),
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
        index_(ctFar),
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

    inline bool farCell() const
    {
        return index_ == ctFar;
    }

    inline int& filterCount()
    {
        return filterCount_;
    }


    inline int filterCount() const
    {
        return filterCount_;
    }


    //- Is the Delaunay cell real, i.e. any real vertex
    inline bool real() const
    {
        return
        (
            this->vertex(0)->real()
         || this->vertex(1)->real()
         || this->vertex(2)->real()
         || this->vertex(3)->real()
        );
    }


    //- Does the Dual vertex form part of a processor patch
    inline bool parallelDualVertex() const
    {
        return
        (
            this->vertex(0)->referred()
         || this->vertex(1)->referred()
         || this->vertex(2)->referred()
         || this->vertex(3)->referred()
        );
    }


    //- Does the Dual vertex form part of a processor patch
    inline Foam::label dualVertexMasterProc() const
    {
        if (!parallelDualVertex())
        {
            return -1;
        }

        // The master processor is the lowest numbered of the four on this tet.

        int masterProc = Foam::Pstream::nProcs() + 1;

        for (int i = 0; i < 4; i++)
        {
            if (this->vertex(i)->referred())
            {
                masterProc = min(masterProc, this->vertex(i)->procIndex());
            }
            else
            {
                masterProc = min(masterProc, Foam::Pstream::myProcNo());
            }
        }

        return masterProc;
    }


    inline Foam::FixedList<Foam::label, 4> processorsAttached() const
    {
        if (!parallelDualVertex())
        {
            return Foam::FixedList<Foam::label, 4>(Foam::Pstream::myProcNo());
        }

        Foam::FixedList<Foam::label, 4> procsAttached
        (
            Foam::Pstream::myProcNo()
        );

        for (int i = 0; i < 4; i++)
        {
            if (this->vertex(i)->referred())
            {
                procsAttached[i] = this->vertex(i)->procIndex();
            }
        }

        return procsAttached;
    }


    // Using the globalIndex object, return a list of four (sorted) global
    // Delaunay vertex indices that uniquely identify this tet in parallel
    inline Foam::FixedList<Foam::label, 4> vertexGlobalIndices
    (
        const Foam::globalIndex& globalDelaunayVertexIndices
    ) const
    {
        // tetVertexGlobalIndices
        Foam::FixedList<Foam::label, 4> tVGI(-1);

        for (int i = 0; i < 4; i++)
        {
            Vertex_handle v = this->vertex(i);

            // Finding the global index of each Delaunay vertex

            if (v->referred())
            {
                // Referred vertices may have negative indices

                tVGI[i] = globalDelaunayVertexIndices.toGlobal
                (
                    v->procIndex(),
                    Foam::mag(v->index())
                );
            }
            else
            {
                tVGI[i] = globalDelaunayVertexIndices.toGlobal
                (
                    Foam::Pstream::myProcNo(),
                    v->index()
                );
            }
        }

        // bubble sort
        for (int i = 0; i < tVGI.size(); i++)
        {
            for (int j = tVGI.size() - 1 ; j > i; j--)
            {
                if (tVGI[j - 1] > tVGI[j])
                {
                    Foam::Swap(tVGI[j - 1], tVGI[j]);
                }
            }
        }

        return tVGI;
    }


    // Is the Delaunay cell part of the final dual mesh, i.e. any vertex form
    // part of the internal or boundary definition
    inline bool internalOrBoundaryDualVertex() const
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
    inline bool anyInternalOrBoundaryDualVertex() const
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


    // Info

        //- Return info proxy.
        //  Used to print indexedCell information to a stream
        Foam::InfoProxy<indexedCell<Gt, Cb> > info() const
        {
            return *this;
        }

        friend Foam::Ostream& Foam::operator<< <Gt, Cb>
        (
            Foam::Ostream&,
            const Foam::InfoProxy<indexedCell<Gt, Cb> >&
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace CGAL

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "indexedCell.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
