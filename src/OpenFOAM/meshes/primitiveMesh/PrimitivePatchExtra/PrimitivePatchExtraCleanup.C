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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Check/fix edges with more than two faces
template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
void Foam::PrimitivePatchExtra<Face, FaceList, PointField, PointType>::
checkEdges
(
    const bool verbose
) const
{
    const labelListList& eFaces = this->edgeFaces();
    const edgeList& edgeLst = this->edges();

    forAll(eFaces, edgeI)
    {
        const labelList& myFaces = eFaces[edgeI];

        // boundary edges have one face
        // interior edges have two faces
        if (myFaces.empty())
        {
            FatalErrorIn
            (
                "PrimitivePatchExtra::checkEdges(bool verbose)"
            )
                << "Edge " << edgeI << " with vertices " << edgeLst[edgeI]
                << " has no edgeFaces"
                << exit(FatalError);
        }
        else if (myFaces.size() > 2)
        {
            WarningIn
            (
                "PrimitivePatchExtra::checkEdges(bool verbose)"
            )
                << "Edge " << edgeI << " with vertices " << edgeLst[edgeI]
                << " has more than 2 faces connected to it : " << myFaces
                << endl;
        }
    }
}


// Check normals and orientation
template
<
    class Face,
    template<class> class FaceList,
    class PointField,
    class PointType
>
Foam::boolList
Foam::PrimitivePatchExtra<Face, FaceList, PointField, PointType>::
checkOrientation
(
    const bool verbose
) const
{
    const FaceList<Face>& faceLst = *this;
    const edgeList& edgeLst = this->edges();
    const labelListList& faceEs = this->faceEdges();
    const label numEdges = this->nEdges();
    const Field<PointType>& pointLst = this->points();
    const vectorField& normLst = this->faceNormals();

    if (ParentType::debug)
    {
        Info<<"checkOrientation:::checkOrientation(bool)" << endl;
    }

    // Check edge normals, face normals, point normals.
    forAll(faceEs, faceI)
    {
        const labelList& edgeLabels = faceEs[faceI];

        if (edgeLabels.size() < 3)
        {
            FatalErrorIn
            (
                "PrimitivePatchExtra::checkOrientation(bool)"
            )
                << "face " << faceLst[faceI]
                << " has fewer than 3 edges. Edges:" << edgeLabels
                << exit(FatalError);
        }

        bool valid = true;
        forAll(edgeLabels, i)
        {
            if (edgeLabels[i] < 0 || edgeLabels[i] >= numEdges)
            {
                WarningIn
                (
                    "PrimitivePatchExtra::checkOrientation(bool)"
                )
                    << "edge number " << edgeLabels[i] << " on face " << faceI
                    << " out-of-range\n"
                    << "This usually means the input surface has "
                    << "edges with more than 2 faces connected.\n"
                    << endl;
                valid = false;
            }
        }
        if (!valid)
        {
            continue;
        }


        //
        //- Compute normal from 3 points, use the first as the origin
        // minor warpage should not be a problem
        const Face& f = faceLst[faceI];
        const point p0(pointLst[f[0]]);
        const point p1(pointLst[f[1]]);
        const point p2(pointLst[f[f.size()-1]]);

        const vector pointNormal((p1 - p0) ^ (p2 - p0));
        if ((pointNormal & normLst[faceI]) < 0)
        {
            FatalErrorIn
            (
                "PrimitivePatchExtra::checkOrientation(bool)"
            )
                << "Normal calculated from points inconsistent with faceNormal"
                << nl
                << "face: " << f << nl
                << "points: " << p0 << ' ' << p1 << ' ' << p2 << nl
                << "pointNormal:" << pointNormal << nl
                << "faceNormal:" << normLst[faceI]
                << exit(FatalError);
        }
    }


    const labelListList& eFaces    = this->edgeFaces();
    const pointField& locPointsLst = this->localPoints();

    // Storage for holding status of edge.
    // True if normal flips across this edge
    boolList borderEdge(numEdges, false);

    forAll(edgeLst, edgeI)
    {
        const edge& e = edgeLst[edgeI];
        const labelList& neighbouringFaces = eFaces[edgeI];

        if (neighbouringFaces.size() == 2)
        {
            // we use localFaces() since edges() are LOCAL
            // these are both already available
            const Face& faceA = this->localFaces()[neighbouringFaces[0]];
            const Face& faceB = this->localFaces()[neighbouringFaces[1]];

            // If the faces are correctly oriented, the edges must go in
            // different directions on connected faces.
            if (faceA.edgeDirection(e) == faceB.edgeDirection(e))
            {
                borderEdge[edgeI] = true;
                if (verbose)
                {
                    WarningIn
                    (
                        "PrimitivePatchExtra::checkOrientation(bool)"
                    )
                        << "face orientation incorrect." << nl
                        << "localEdge[" << edgeI << "] " << e
                        << " between faces:" << nl
                        << "  face[" << neighbouringFaces[0] << "] "
                        << faceLst[neighbouringFaces[0]]
                        << "   localFace: " << faceA
                        << nl
                        << "  face[" << neighbouringFaces[1] << "] "
                        << faceLst[neighbouringFaces[1]]
                        << "   localFace: " << faceB
                        << endl;
                }
            }
        }
        else if (neighbouringFaces.size() != 1)
        {
            if (verbose)
            {
                WarningIn
                (
                    "PrimitivePatchExtra::checkOrientation(bool)"
                )
                    << "Wrong number of edge neighbours." << nl
                    << "edge[" << edgeI << "] " << e
                    << " with points:" << locPointsLst[e.start()]
                    << ' ' << locPointsLst[e.end()]
                    << " has neighbouringFaces:" << neighbouringFaces << endl;
            }
            borderEdge[edgeI] = true;
        }
    }

    return borderEdge;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
