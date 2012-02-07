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

template<class Gt, class Vb>
inline CGAL::indexedVertex<Gt, Vb>::indexedVertex()
:
    Vb(),
    index_(vtInternal),
    type_(vtInternal),
    alignment_(),
    targetCellSize_(0.0),
    vertexFixed_(false)
{}


template<class Gt, class Vb>
inline CGAL::indexedVertex<Gt, Vb>::indexedVertex(const Point& p)
:
    Vb(p),
    index_(vtInternal),
    type_(vtInternal),
    alignment_(),
    targetCellSize_(0.0),
    vertexFixed_(false)
{}


template<class Gt, class Vb>
inline CGAL::indexedVertex<Gt, Vb>::indexedVertex
(
    const Point& p,
    int index,
    int type
)
:
    Vb(p),
    index_(index),
    type_(type),
    alignment_(),
    targetCellSize_(0.0),
    vertexFixed_(false)
{}


template<class Gt, class Vb>
inline CGAL::indexedVertex<Gt, Vb>::indexedVertex(const Point& p, Cell_handle f)
:
    Vb(f, p),
    index_(vtInternal),
    type_(vtInternal),
    alignment_(),
    targetCellSize_(0.0),
    vertexFixed_(false)
{}


template<class Gt, class Vb>
inline CGAL::indexedVertex<Gt, Vb>::indexedVertex(Cell_handle f)
:
    Vb(f),
    index_(vtInternal),
    type_(vtInternal),
    alignment_(),
    targetCellSize_(0.0),
    vertexFixed_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Gt, class Vb>
inline int& CGAL::indexedVertex<Gt, Vb>::index()
{
    return index_;
}


template<class Gt, class Vb>
inline int CGAL::indexedVertex<Gt, Vb>::index() const
{
    return index_;
}


template<class Gt, class Vb>
inline int& CGAL::indexedVertex<Gt, Vb>::type()
{
    return type_;
}


template<class Gt, class Vb>
inline int CGAL::indexedVertex<Gt, Vb>::type() const
{
    return type_;
}


template<class Gt, class Vb>
inline Foam::tensor& CGAL::indexedVertex<Gt, Vb>::alignment()
{
    return alignment_;
}


template<class Gt, class Vb>
inline const Foam::tensor& CGAL::indexedVertex<Gt, Vb>::alignment() const
{
    return alignment_;
}


template<class Gt, class Vb>
inline Foam::scalar& CGAL::indexedVertex<Gt, Vb>::targetCellSize()
{
    return targetCellSize_;
}


template<class Gt, class Vb>
inline Foam::scalar CGAL::indexedVertex<Gt, Vb>::targetCellSize() const
{
    return targetCellSize_;
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::uninitialised() const
{
    return type_ == vtInternal && index_ == vtInternal;
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::farPoint() const
{
    return type_ == vtFar;
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::internalPoint() const
{
    return internalPoint(type_);
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::internalPoint(int type)
{
    return type <= vtInternal;
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::referred() const
{
    return (type_ < 0 && type_ > vtFar);
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::referredInternalOrBoundaryPoint() const
{
    return referred() && index_ >= 0;
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::referredExternal() const
{
    return referred() && index_ < 0;
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::real() const
{
    return internalPoint() || pairPoint();
}


template<class Gt, class Vb>
inline int CGAL::indexedVertex<Gt, Vb>::procIndex() const
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


template<class Gt, class Vb>
inline int CGAL::indexedVertex<Gt, Vb>::encodeProcIndex(int procI)
{
    return -(procI + 1);
}


template<class Gt, class Vb>
inline void CGAL::indexedVertex<Gt, Vb>::setInternal()
{
    type_ = vtInternal;
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::nearBoundary() const
{
    return type_ == vtNearBoundary;
}


template<class Gt, class Vb>
inline void CGAL::indexedVertex<Gt, Vb>::setNearBoundary()
{
    type_ = vtNearBoundary;
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::pairPoint() const
{
    return type_ >= 0;
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::ppMaster() const
{
    return ppMaster(index_, type_);
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::ppMaster(int index, int type)
{
    if (index >= 0 && type > index)
    {
        return true;
    }

    return false;
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::ppSlave() const
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


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::internalOrBoundaryPoint() const
{
    return internalOrBoundaryPoint(index_, type_);
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::internalOrBoundaryPoint
(
    int index,
    int type
)
{
    return internalPoint(type) || ppMaster(index, type);
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::nearOrOnBoundary() const
{
    return pairPoint() || nearBoundary();
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::anyInternalOrBoundaryPoint() const
{
    return internalOrBoundaryPoint() || referredInternalOrBoundaryPoint();
}


template<class Gt, class Vb>
inline bool CGAL::indexedVertex<Gt, Vb>::isVertexFixed() const
{
    return vertexFixed_;
}


template<class Gt, class Vb>
inline void CGAL::indexedVertex<Gt, Vb>::setVertexFixed() const
{
    vertexFixed_ = true;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //
