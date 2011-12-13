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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Gt, class Cb>
CGAL::indexedCell<Gt, Cb>::indexedCell()
:
    Cb(),
    index_(ctFar),
    filterCount_(0)
{}


template<class Gt, class Cb>
CGAL::indexedCell<Gt, Cb>::indexedCell
(
    Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Vertex_handle v3
)
:
    Cb(v0, v1, v2, v3),
    index_(ctFar),
    filterCount_(0)
{}


template<class Gt, class Cb>
CGAL::indexedCell<Gt, Cb>::indexedCell
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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Gt, class Cb>
int& CGAL::indexedCell<Gt, Cb>::cellIndex()
{
    return index_;
}


template<class Gt, class Cb>
int CGAL::indexedCell<Gt, Cb>::cellIndex() const
{
    return index_;
}


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::farCell() const
{
    return index_ == ctFar;
}


template<class Gt, class Cb>
inline int& CGAL::indexedCell<Gt, Cb>::filterCount()
{
    return filterCount_;
}


template<class Gt, class Cb>
inline int CGAL::indexedCell<Gt, Cb>::filterCount() const
{
    return filterCount_;
}


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::real() const
{
    return
    (
        this->vertex(0)->real()
     || this->vertex(1)->real()
     || this->vertex(2)->real()
     || this->vertex(3)->real()
    );
}


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::parallelDualVertex() const
{
    return
    (
        this->vertex(0)->referred()
     || this->vertex(1)->referred()
     || this->vertex(2)->referred()
     || this->vertex(3)->referred()
    );
}


template<class Gt, class Cb>
inline Foam::label CGAL::indexedCell<Gt, Cb>::dualVertexMasterProc() const
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


template<class Gt, class Cb>
inline Foam::FixedList<Foam::label, 4>
CGAL::indexedCell<Gt, Cb>::processorsAttached() const
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


template<class Gt, class Cb>
inline Foam::tetCell CGAL::indexedCell<Gt, Cb>::vertexGlobalIndices
(
    const Foam::globalIndex& globalDelaunayVertexIndices
) const
{
    // tetVertexGlobalIndices
    Foam::tetCell tVGI;

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


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::internalOrBoundaryDualVertex() const
{
    return
    (
        this->vertex(0)->internalOrBoundaryPoint()
     || this->vertex(1)->internalOrBoundaryPoint()
     || this->vertex(2)->internalOrBoundaryPoint()
     || this->vertex(3)->internalOrBoundaryPoint()
    );
}


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::anyInternalOrBoundaryDualVertex() const
{
    return
    (
        this->vertex(0)->anyInternalOrBoundaryPoint()
     || this->vertex(1)->anyInternalOrBoundaryPoint()
     || this->vertex(2)->anyInternalOrBoundaryPoint()
     || this->vertex(3)->anyInternalOrBoundaryPoint()
    );
}


template<class Gt, class Cb>
inline bool CGAL::indexedCell<Gt, Cb>::boundaryDualVertex() const
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


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //
