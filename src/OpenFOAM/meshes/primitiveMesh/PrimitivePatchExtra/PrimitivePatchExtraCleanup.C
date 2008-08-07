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

\*---------------------------------------------------------------------------*/

#include "PrimitivePatchExtra.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Check/fix edges with more than two faces
template
<
    class Face,
    template<class> class ListType,
    class PointField,
    class PointType
>
void Foam::PrimitivePatchExtra<Face, ListType, PointField, PointType>::
checkEdges
(
    const bool verbose
) const
{
    const labelListList& eFaces =
       PrimitivePatch<Face, ListType, PointField>::edgeFaces();

    const edgeList& edgeLst =
        PrimitivePatch<Face, ListType, PointField>::edges();

    forAll (eFaces, edgeI)
    {
        const labelList& myFaces = eFaces[edgeI];

        if (myFaces.size() == 0)
        {
            FatalErrorIn("PrimitivePatchExtra::checkEdges(bool verbose)")
                << "Edge " << edgeI << " with vertices " << edgeLst[edgeI]
                << " has no edgeFaces"
                << exit(FatalError);
        }
        else if (myFaces.size() > 2)
        {
            WarningIn
            (
                "PrimitivePatchExtra::checkEdges(bool verbose)"
            )   << "Edge " << edgeI << " with vertices " << edgeLst[edgeI]
                << " has more than 2 faces connected to it : " << myFaces
                << endl;
        }
    }
}


// Check normals and orientation
template
<
    class Face,
    template<class> class ListType,
    class PointField,
    class PointType
>
Foam::boolList
Foam::PrimitivePatchExtra<Face, ListType, PointField, PointType>::
checkOrientation
(
    const bool verbose
) const
{
    const edgeList& edgeLst = TemplateType::edges();
    const labelListList& faceEs = TemplateType::faceEdges();
    const List<FaceType>& faceLst = TemplateType::faces();
    const label numEdges = TemplateType::nEdges();
    const pointField& pointLst = TemplateType::points();
    const vectorField& normLst = TemplateType::faceNormals();

    // Check edge normals, face normals, point normals.
    forAll (faceEs, faceI)
    {
        const labelList& edgeLabels = faceEs[faceI];

        if (edgeLabels.size() < 3)
        {
            FatalErrorIn("PrimitivePatchExtra::checkOrientation(bool)")
                << "face " << faceLst[faceI]
                << " has fewer than 3 edges. Edges:" << edgeLabels
                << exit(FatalError);
        }

        bool valid = true;
        forAll (edgeLabels, i)
        {
            if (edgeLabels[i] < 0 || edgeLabels[i] >= numEdges)
            {
                WarningIn
                (
                    "PrimitivePatchExtra::checkOrientation(bool)"
                )   << "edge number " << edgeLabels[i] << " on face " << faceI
                    << " out of range"
                    << "\nThis usually means that the input surface has "
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
        //
        const FaceType& f = faceLst[faceI];
        const point p0(pointLst[f[0]]);
        const point p1(pointLst[f[1]]);
        const point p2(pointLst[f[f.size()-1]]);

        const vector pointNormal((p1 - p0) ^ (p2 - p0));
        if ((pointNormal & normLst[faceI]) < 0)
        {
            FatalErrorIn("PrimitivePatchExtra::checkOrientation(bool)")
                << "Normal calculated from points not consistent with"
                " faceNormal" << endl
                << "face: " << f << endl
                << "points: " << p0 << ' ' << p1 << ' ' << p2 << endl
                << "pointNormal:" << pointNormal << endl
                << "faceNormal:" << normLst[faceI]
                << exit(FatalError);
        }
    }


    const labelListList& eFaces = TemplateType::edgeFaces();
    const pointField& locPointsLst = TemplateType::localPoints();

    // Storage for holding status of edge. True if normal flips across this
    // edge
    boolList borderEdge(numEdges, false);

    forAll (edgeLst, edgeI)
    {
        const labelList& neighbours = eFaces[edgeI];

        if (neighbours.size() == 2)
        {
            // Two faces, A and B. Check if edge orientation is
            // anticlockwise on both.
            const labelList& fEdgesA = faceEs[neighbours[0]];
            const labelList& fEdgesB = faceEs[neighbours[1]];

            // Get next edge after edgeI
            label nextEdgeA = fEdgesA.fcIndex(findIndex(fEdgesA, edgeI));
            label nextEdgeB = fEdgesB.fcIndex(findIndex(fEdgesB, edgeI));

            // Now check if nextEdgeA and nextEdgeB have any common points
            if
            (
                edgeLst[nextEdgeA].start() == edgeLst[nextEdgeB].start()
             || edgeLst[nextEdgeA].start() == edgeLst[nextEdgeB].end()
             || edgeLst[nextEdgeA].end() == edgeLst[nextEdgeB].start()
             || edgeLst[nextEdgeA].end() == edgeLst[nextEdgeB].end()
            )
            {
                borderEdge[edgeI] = true;
                if (verbose)
                {
                    // just list first three points
                    // to simplify generating the message
                    WarningIn("PrimitivePatchExtra::checkOrientation(bool)")
                        << "face orientation incorrect." << nl
                        << "edge neighbours:" << neighbours << nl
                        << "face " << neighbours[0] << " has edges "
                        << fEdgesA << nl
                        << "    with points " << nl
                        << "    " << edgeLst[fEdgesA[0]].start() << ' '
                        << edgeLst[fEdgesA[0]].end() << nl
                        << "    " << edgeLst[fEdgesA[1]].start() << ' '
                        << edgeLst[fEdgesA[1]].end() << nl
                        << "    " << edgeLst[fEdgesA[2]].start() << ' '
                        << edgeLst[fEdgesA[2]].end()
                        << endl
                        << "face " << neighbours[1] << " has edges "
                        << fEdgesB << nl
                        << "    with points " << nl
                        << "    " << edgeLst[fEdgesB[0]].start() << ' '
                        << edgeLst[fEdgesB[0]].end() << nl
                        << "    " << edgeLst[fEdgesB[1]].start() << ' '
                        << edgeLst[fEdgesB[1]].end() << nl
                        << "    " << edgeLst[fEdgesB[2]].start() << ' '
                        << edgeLst[fEdgesB[2]].end() << nl
                        << endl;
                }
            }
        }
        else if (neighbours.size() != 1)
        {
            if (verbose)
            {
                const edge& e = edgeLst[edgeI];
                WarningIn("PrimitivePatchExtra::checkOrientation(bool)")
                    << "Wrong number of edge neighbours." << endl
                    << "Edge:" << e
                    << "with points:" << locPointsLst[e.start()]
                    << ' ' << locPointsLst[e.end()]
                    << " has neighbours:" << neighbours << endl;
            }
            borderEdge[edgeI] = true;
        }
    }

    return borderEdge;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
