/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "PrimitivePatchExtra.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
Foam::PrimitivePatchExtra<Face, FaceList, PointField, PointType>::
PrimitivePatchExtra
(
    const FaceList<Face>& faces,
    const Field<PointType>& points
)
:
    ParentType(faces, points),
    sortedEdgeFacesPtr_(NULL),
    edgeOwnerPtr_(NULL)
{}


// Construct as copy
template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
Foam::PrimitivePatchExtra<Face, FaceList, PointField, PointType>::
PrimitivePatchExtra
(
    const PrimitivePatchExtra<Face, FaceList, PointField, PointType>& pp
)
:
    ParentType(pp),
    sortedEdgeFacesPtr_(NULL),
    edgeOwnerPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
Foam::PrimitivePatchExtra<Face, FaceList, PointField, PointType>::
~PrimitivePatchExtra()
{
    clearOut();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
void Foam::PrimitivePatchExtra<Face, FaceList, PointField, PointType>::
clearOut()
{
    ParentType::clearOut();
    clearTopology();
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
void Foam::PrimitivePatchExtra<Face, FaceList, PointField, PointType>::
clearTopology()
{
    ParentType::clearTopology();
    deleteDemandDrivenData(sortedEdgeFacesPtr_);
    deleteDemandDrivenData(edgeOwnerPtr_);
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
const Foam::labelListList&
Foam::PrimitivePatchExtra<Face, FaceList, PointField, PointType>::
sortedEdgeFaces() const
{
    if (!sortedEdgeFacesPtr_)
    {
        calcSortedEdgeFaces();
    }

    return *sortedEdgeFacesPtr_;
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
const Foam::labelList&
Foam::PrimitivePatchExtra<Face, FaceList, PointField, PointType>::
edgeOwner() const
{
    if (!edgeOwnerPtr_)
    {
        calcEdgeOwner();
    }

    return *edgeOwnerPtr_;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PrimitivePatchExtraAddressing.C"
#include "PrimitivePatchExtraCleanup.C"
#include "PrimitivePatchExtraSearch.C"

// ************************************************************************* //
