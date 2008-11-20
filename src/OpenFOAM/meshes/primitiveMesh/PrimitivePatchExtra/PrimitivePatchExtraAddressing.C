/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

Description
    Contains fix for PrimitivePatch addressing (which doesn't work if surface
    is non-manifold). Should be moved into PrimitivePatch.

\*---------------------------------------------------------------------------*/

#include "PrimitivePatchExtra.H"
#include "HashTable.H"
#include "SortableList.H"
#include "transform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
void Foam::PrimitivePatchExtra<Face, FaceList, PointField, PointType>::
calcSortedEdgeFaces() const
{
    if (sortedEdgeFacesPtr_)
    {
        FatalErrorIn
        (
            "PrimitivePatchExtra<Face, FaceList, PointField>::"
            "calcSortedEdgeFaces()"
        )
            << "sortedEdgeFacesPtr_ already set"
            << abort(FatalError);
    }

    const labelListList& eFaces = this->edgeFaces();
    const edgeList& edgeLst = this->edges();
    const Field<PointType>& locPointLst = this->localPoints();
    const List<Face>& locFaceLst  = this->localFaces();

    // create the lists for the various results. (resized on completion)
    sortedEdgeFacesPtr_ = new labelListList(eFaces.size());
    labelListList& sortedEdgeFaces = *sortedEdgeFacesPtr_;

    forAll(eFaces, edgeI)
    {
        const labelList& myFaceNbs = eFaces[edgeI];

        if (myFaceNbs.size() > 2)
        {
            // Get point on edge and normalized direction of edge (= e2 base
            // of our coordinate system)
            const edge& e = edgeLst[edgeI];

            const point& edgePt = locPointLst[e.start()];

            vector e2 = e.vec(locPointLst);
            e2 /= mag(e2) + VSMALL;

            // Get opposite vertex for 0th face
            const Face& f = locFaceLst[myFaceNbs[0]];

            label fp0 = findIndex(f, e[0]);
            label fp1 = f.fcIndex(fp0);
            label vertI = (f[fp1] != e[1] ? f[fp1] : f.fcIndex(fp1));

            // Get vector normal both to e2 and to edge from opposite vertex
            // to edge (will be x-axis of our coordinate system)
            vector e0 = e2 ^ (locPointLst[vertI] - edgePt);
            e0 /= mag(e0) + VSMALL;

            // Get y-axis of coordinate system
            vector e1 = e2 ^ e0;

            SortableList<scalar> faceAngles(myFaceNbs.size());

            // e0 is reference so angle is 0
            faceAngles[0] = 0;

            for (label nbI = 1; nbI < myFaceNbs.size(); nbI++)
            {
                // Get opposite vertex
                const Face& f = locFaceLst[myFaceNbs[nbI]];
                label fp0 = findIndex(f, e[0]);
                label fp1 = f.fcIndex(fp0);
                label vertI = (f[fp1] != e[1] ? f[fp1] : f.fcIndex(fp1));

                vector vec = e2 ^ (locPointLst[vertI] - edgePt);
                vec /= mag(vec) + VSMALL;

                faceAngles[nbI] = pseudoAngle
                (
                    e0,
                    e1,
                    vec
                );
            }

            faceAngles.sort();

            sortedEdgeFaces[edgeI] = IndirectList<label>
            (
                myFaceNbs,
                faceAngles.indices()
            );
        }
        else
        {
            // No need to sort. Just copy.
            sortedEdgeFaces[edgeI] = myFaceNbs;
        }
    }
}


template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
void Foam::PrimitivePatchExtra<Face, FaceList, PointField, PointType>::
calcEdgeOwner() const
{
    if (edgeOwnerPtr_)
    {
        FatalErrorIn
        (
            "PrimitivePatchExtra<Face, FaceList, PointField>::"
            "calcEdgeOwner()"
        )
            << "edgeOwnerPtr_ already set"
            << abort(FatalError);
    }

    // create the owner list
    edgeOwnerPtr_ = new labelList(this->nEdges());
    labelList& edgeOwner = *edgeOwnerPtr_;

    const edgeList& edgeLst = this->edges();
    const labelListList& eFaces = this->edgeFaces();
    const List<Face>& locFaceLst = this->localFaces();


    forAll(edgeLst, edgeI)
    {
        const edge& e = edgeLst[edgeI];

        const labelList& myFaces = eFaces[edgeI];

        if (myFaces.size() == 1)
        {
            edgeOwner[edgeI] = myFaces[0];
        }
        else
        {
            // Find the first face whose vertices are aligned with the edge.
            // (in case of multiply connected edge the best we can do)
            edgeOwner[edgeI] = -1;

            forAll(myFaces, i)
            {
                const Face& f = locFaceLst[myFaces[i]];

                if (f.findEdge(e) > 0)
                {
                    edgeOwner[edgeI] = myFaces[i];
                    break;
                }
            }

            if (edgeOwner[edgeI] == -1)
            {
                FatalErrorIn
                (
                    "PrimitivePatchExtra<Face, FaceList, PointField>::"
                    "calcEdgeOwner()"
                )
                    << "Edge " << edgeI << " vertices:" << e
                    << " is used by faces " << myFaces
                    << " vertices:"
                    << IndirectList<Face>(locFaceLst, myFaces)()
                    << " none of which use the edge vertices in the same order"
                    << nl << "I give up" << abort(FatalError);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
