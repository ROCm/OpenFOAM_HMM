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

\*---------------------------------------------------------------------------*/

//#include "PrimitivePatch.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Face, template<class> class FaceList, class PointField>
PrimitivePatch<Face, FaceList, PointField>::PrimitivePatch
(
    const FaceList<Face>& faces,
    const pointField& points
)
:
    FaceList<Face>(faces),
    points_(points),
    edgesPtr_(NULL),
    nInternalEdges_(-1),
    boundaryPointsPtr_(NULL),
    faceFacesPtr_(NULL),
    edgeFacesPtr_(NULL),
    faceEdgesPtr_(NULL),
    pointEdgesPtr_(NULL),
    pointFacesPtr_(NULL),
    localFacesPtr_(NULL),
    meshPointsPtr_(NULL),
    meshPointMapPtr_(NULL),
    edgeLoopsPtr_(NULL),
    localPointsPtr_(NULL),
    localPointOrderPtr_(NULL),
    faceNormalsPtr_(NULL),
    pointNormalsPtr_(NULL)
{}


// Construct as copy
template<class Face, template<class> class FaceList, class PointField>
PrimitivePatch<Face, FaceList, PointField>::PrimitivePatch
(
    const PrimitivePatch<Face, FaceList, PointField>& pp
)
:
    PrimitivePatchName(),
    FaceList<Face>(pp),
    points_(pp.points_),
    edgesPtr_(NULL),
    nInternalEdges_(-1),
    boundaryPointsPtr_(NULL),
    faceFacesPtr_(NULL),
    edgeFacesPtr_(NULL),
    faceEdgesPtr_(NULL),
    pointEdgesPtr_(NULL),
    pointFacesPtr_(NULL),
    localFacesPtr_(NULL),
    meshPointsPtr_(NULL),
    meshPointMapPtr_(NULL),
    edgeLoopsPtr_(NULL),
    localPointsPtr_(NULL),
    localPointOrderPtr_(NULL),
    faceNormalsPtr_(NULL),
    pointNormalsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Face, template<class> class FaceList, class PointField>
PrimitivePatch<Face, FaceList, PointField>::~PrimitivePatch()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Correct patch after moving points
template<class Face, template<class> class FaceList, class PointField>
void PrimitivePatch<Face, FaceList, PointField>::movePoints(const pointField&)
{
    if (debug)
    {
        Info<< "PrimitivePatch<Face, FaceList, PointField>::movePoints() : "
            << "recalculating PrimitivePatch geometry following mesh motion"
            << endl;
    }

    clearGeom();
}


template<class Face, template<class> class FaceList, class PointField>
const edgeList&
PrimitivePatch<Face, FaceList, PointField>::edges() const
{
    if (!edgesPtr_)
    {
        calcAddressing();
    }

    return *edgesPtr_;
}


template<class Face, template<class> class FaceList, class PointField>
label PrimitivePatch<Face, FaceList, PointField>::nInternalEdges() const
{
    if (!edgesPtr_)
    {
        calcAddressing();
    }

    return nInternalEdges_;
}


template<class Face, template<class> class FaceList, class PointField>
const labelList&
PrimitivePatch<Face, FaceList, PointField>::boundaryPoints() const
{
    if (!boundaryPointsPtr_)
    {
        calcBdryPoints();
    }

    return *boundaryPointsPtr_;
}


template<class Face, template<class> class FaceList, class PointField>
const labelListList&
PrimitivePatch<Face, FaceList, PointField>::faceFaces() const
{
    if (!faceFacesPtr_)
    {
        calcAddressing();
    }

    return *faceFacesPtr_;
}


template<class Face, template<class> class FaceList, class PointField>
const labelListList&
PrimitivePatch<Face, FaceList, PointField>::edgeFaces() const
{
    if (!edgeFacesPtr_)
    {
        calcAddressing();
    }

    return *edgeFacesPtr_;
}


template<class Face, template<class> class FaceList, class PointField>
const labelListList&
PrimitivePatch<Face, FaceList, PointField>::faceEdges() const
{
    if (!faceEdgesPtr_)
    {
        calcAddressing();
    }

    return *faceEdgesPtr_;
}


template<class Face, template<class> class FaceList, class PointField>
const labelListList&
PrimitivePatch<Face, FaceList, PointField>::pointEdges() const
{
    if (!pointEdgesPtr_)
    {
        calcPointEdges();
    }

    return *pointEdgesPtr_;
}


template<class Face, template<class> class FaceList, class PointField>
const labelListList&
PrimitivePatch<Face, FaceList, PointField>::pointFaces() const
{
    if (!pointFacesPtr_)
    {
        calcPointFaces();
    }

    return *pointFacesPtr_;
}


template<class Face, template<class> class FaceList, class PointField>
const List<Face>&
PrimitivePatch<Face, FaceList, PointField>::localFaces() const
{
    if (!localFacesPtr_)
    {
        calcMeshData();
    }

    return *localFacesPtr_;
}


template<class Face, template<class> class FaceList, class PointField>
const labelList&
PrimitivePatch<Face, FaceList, PointField>::meshPoints() const
{
    if (!meshPointsPtr_)
    {
        calcMeshData();
    }

    return *meshPointsPtr_;
}


template<class Face, template<class> class FaceList, class PointField>
const Map<label>&
PrimitivePatch<Face, FaceList, PointField>::meshPointMap() const
{
    if (!meshPointMapPtr_)
    {
        calcMeshPointMap();
    }

    return *meshPointMapPtr_;
}


template<class Face, template<class> class FaceList, class PointField>
const pointField&
PrimitivePatch<Face, FaceList, PointField>::localPoints() const
{
    if (!localPointsPtr_)
    {
        calcLocalPoints();
    }

    return *localPointsPtr_;
}


template<class Face, template<class> class FaceList, class PointField>
const labelList&
PrimitivePatch<Face, FaceList, PointField>::localPointOrder() const
{
    if (!localPointOrderPtr_)
    {
        calcLocalPointOrder();
    }

    return *localPointOrderPtr_;
}


template<class Face, template<class> class FaceList, class PointField>
label PrimitivePatch<Face, FaceList, PointField>::whichPoint
(
    const label gp
) const
{
    Map<label>::const_iterator gpIter = meshPointMap().find(gp);

    if (gpIter != meshPointMap().end())
    {
        return gpIter();
    }
    else
    {
        // Not found
        return -1;
    }
}


template<class Face, template<class> class FaceList, class PointField>
const vectorField&
PrimitivePatch<Face, FaceList, PointField>::faceNormals() const
{
    if (!faceNormalsPtr_)
    {
        calcFaceNormals();
    }

    return *faceNormalsPtr_;
}


template<class Face, template<class> class FaceList, class PointField>
const vectorField&
PrimitivePatch<Face, FaceList, PointField>::pointNormals() const
{
    if (!pointNormalsPtr_)
    {
        calcPointNormals();
    }

    return *pointNormalsPtr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Face, template<class> class FaceList, class PointField>
void PrimitivePatch<Face, FaceList, PointField>::operator=
(
    const PrimitivePatch<Face, FaceList, PointField>& pp
)
{
    clearOut();

    FaceList<Face>::operator=(pp);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PrimitivePatchAddressing.C"
#include "PrimitivePatchEdgeLoops.C"
#include "PrimitivePatchClear.C"
#include "PrimitivePatchBdryPoints.C"
#include "PrimitivePatchLocalPointOrder.C"
#include "PrimitivePatchMeshData.C"
#include "PrimitivePatchMeshEdges.C"
#include "PrimitivePatchPointAddressing.C"
#include "PrimitivePatchProjectPoints.C"
#include "PrimitivePatchCheck.C"

// ************************************************************************* //
