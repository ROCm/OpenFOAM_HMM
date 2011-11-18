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
    under the terms_ of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "shortEdgeFilter2D.H"

namespace Foam
{

defineTypeNameAndDebug(shortEdgeFilter2D, 0);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shortEdgeFilter2D::shortEdgeFilter2D
(
    const Foam::CV2D& cv2Dmesh,
    const dictionary& dict
)
:
    cv2Dmesh_(cv2Dmesh),
    shortEdgeFilterFactor_(readScalar(dict.lookup("shortEdgeFilterFactor"))),
    edgeAttachedToBoundaryFactor_
    (
        dict.lookupOrDefault<scalar>("edgeAttachedToBoundaryFactor", 2.0)
    ),
    patchNames_(wordList()),
    patchSizes_(labelList()),
    mapEdgesRegion_()
{
    point2DField points2D;
    faceList faces;

    cv2Dmesh.calcDual
    (
        points2D,
        faces,
        patchNames_,
        patchSizes_,
        mapEdgesRegion_
    );

    pointField points(points2D.size());
    forAll(points, ip)
    {
        points[ip] = cv2Dmesh.toPoint3D(points2D[ip]);
    }

    points2D.clear();

    ms_ = MeshedSurface<face>(xferMove(points), xferMove(faces));

    Info<< "Meshed surface stats before edge filtering :" << endl;
    ms_.writeStats(Info);

    if (debug)
    {
        writeInfo(Info);

        ms_.write("MeshedSurface_preFilter.obj");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::shortEdgeFilter2D::~shortEdgeFilter2D()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void
Foam::shortEdgeFilter2D::filter()
{
    // These are global indices.
    const pointField& points = ms_.points();
    const edgeList& edges = ms_.edges();
    const faceList& faces = ms_.faces();
    const labelList& meshPoints = ms_.meshPoints();
    const labelList& boundaryPoints = ms_.boundaryPoints();

    label maxChain = 0;
    label nPointsToRemove = 0;

    labelList pointsToRemove(ms_.points().size(), -1);

    // List of number of vertices in a face.
    labelList newFaceVertexCount(faces.size(), -1);
    forAll(faces, faceI)
    {
        newFaceVertexCount[faceI] = faces[faceI].size();
    }

    // Check if the point is a boundary point. Flag if it is so that
    // it will not be deleted.
    boolList boundaryPointFlags(points.size(), false);
    // This has been removed, otherwise small edges on the boundary are not
    // removed.
    /*  forAll(boundaryPointFlags, pointI)
    {
        forAll(boundaryPoints, bPoint)
        {
            if (meshPoints[boundaryPoints[bPoint]] == pointI)
            {
                boundaryPointFlags[pointI] = true;
            }
        }
    }*/

    // Check if an edge has a boundary point. It it does the edge length
    // will be doubled when working out its length.
    Info<< "    Marking edges attached to boundaries." << endl;
    boolList edgeAttachedToBoundary(edges.size(), false);
    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];
        const label startVertex = e.start();
        const label endVertex = e.end();

        forAll(boundaryPoints, bPoint)
        {
            if
            (
                boundaryPoints[bPoint] == startVertex
             || boundaryPoints[bPoint] == endVertex
            )
            {
                edgeAttachedToBoundary[edgeI] = true;
            }
        }
    }

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        // get the vertices of that edge.
        const label startVertex = e.start();
        const label endVertex = e.end();

        scalar edgeLength =
            mag
            (
                points[meshPoints[e.start()]]
              - points[meshPoints[e.end()]]
            );

        if (edgeAttachedToBoundary[edgeI])
        {
            edgeLength *= edgeAttachedToBoundaryFactor_;
        }

        scalar shortEdgeFilterValue = 0.0;

        const labelList& psEdges = ms_.pointEdges()[startVertex];
        const labelList& peEdges = ms_.pointEdges()[endVertex];

        forAll(psEdges, psEdgeI)
        {
            const edge& psE = edges[psEdges[psEdgeI]];
            if (edgeI != psEdges[psEdgeI])
            {
                shortEdgeFilterValue +=
                    mag
                    (
                        points[meshPoints[psE.start()]]
                       -points[meshPoints[psE.end()]]
                    );
            }
        }

        forAll(peEdges, peEdgeI)
        {
            const edge& peE = edges[peEdges[peEdgeI]];
            if (edgeI != peEdges[peEdgeI])
            {
                shortEdgeFilterValue +=
                    mag
                    (
                        points[meshPoints[peE.start()]]
                       -points[meshPoints[peE.end()]]
                    );
            }
        }

        shortEdgeFilterValue *=
            shortEdgeFilterFactor_
           /(psEdges.size() + peEdges.size() - 2);

        if (edgeLength < shortEdgeFilterValue)
        {
            bool flagDegenerateFace = false;
            const labelList& pFaces = ms_.pointFaces()[startVertex];

            forAll(pFaces, pFaceI)
            {
                const face& f = ms_.localFaces()[pFaces[pFaceI]];
                forAll(f, fp)
                {
                    // If the edge is part of this face...
                    if (f[fp] == endVertex)
                    {
                        // If deleting vertex would create a triangle, don't!
                        if (newFaceVertexCount[pFaces[pFaceI]] < 4)
                        {
                            flagDegenerateFace = true;
                        }
                        else
                        {
                            newFaceVertexCount[pFaces[pFaceI]]--;
                        }
                    }
                    // If the edge is not part of this face...
                    else
                    {
                        // Deleting vertex of a triangle...
                        if (newFaceVertexCount[pFaces[pFaceI]] < 3)
                        {
                            flagDegenerateFace = true;
                        }
                    }
                }
            }

            // This if statement determines whether a point should be deleted.
            if
            (
                pointsToRemove[meshPoints[startVertex]] == -1
             && pointsToRemove[meshPoints[endVertex]] == -1
             && !boundaryPointFlags[meshPoints[startVertex]]
             && !flagDegenerateFace
            )
            {
                pointsToRemove[meshPoints[startVertex]] =
                    meshPoints[endVertex];
                ++nPointsToRemove;
            }
        }
    }

    label totalNewPoints = points.size() - nPointsToRemove;

    pointField newPoints(totalNewPoints, vector(0, 0, 0));
    labelList newPointNumbers(points.size(), -1);
    label numberRemoved=0;

    forAll(points, pointI)
    {
        // If the point is NOT going to be removed.
        if (pointsToRemove[pointI] == -1)
        {
            newPoints[pointI-numberRemoved] = points[pointI];
            newPointNumbers[pointI] =  pointI-numberRemoved;
        }
        else
        {
            numberRemoved++;
        }
    }

    // Need a new faceList
    faceList newFaces(faces.size());
    label newFaceI = 0;

    labelList newFace;
    label newFaceSize = 0;

    // Now need to iterate over the faces and remove points. Global index.
    forAll(faces, faceI)
    {
        const face& f = faces[faceI];

        newFace.clear();
        newFace.setSize(f.size());
        newFaceSize = 0;

        forAll(f, fp)
        {
            label pointI = f[fp];
            // If not removing the point, then add it to the new face.
            if (pointsToRemove[pointI] == -1)
            {
                newFace[newFaceSize++] = newPointNumbers[pointI];
            }
            else
            {
                label newPointI = pointsToRemove[pointI];
                // Replace deleted point with point that it is being
                // collapsed to.
                if
                (
                    f.nextLabel(fp) != newPointI
                 && f.prevLabel(fp) != newPointI
                )
                {
                    label pChain = newPointI;
                    label totalChain = 0;
                    for (label nChain = 0; nChain <= totalChain; ++nChain)
                    {
                        if (newPointNumbers[pChain] != -1)
                        {
                            newFace[newFaceSize++] = newPointNumbers[pChain];
                            newPointNumbers[pointI]
                                = newPointNumbers[pChain];
                            maxChain = max(totalChain, maxChain);
                        }
                        else
                        {
                            WarningIn("shortEdgeFilter")
                                << "Point " << pChain
                                << " marked for deletion as well as point "
                                << pointI << nl
                                << "    Incrementing maxChain by 1 from "
                                << totalChain << " to " << totalChain + 1
                                << endl;
                            totalChain++;
                        }
                        pChain = pointsToRemove[pChain];
                    }
                }
                else
                {
                    if (newPointNumbers[newPointI] != -1)
                    {
                        newPointNumbers[pointI] = newPointNumbers[newPointI];
                    }
                }
            }
        }

        newFace.setSize(newFaceSize);

        if (newFace.size() > 2)
        {
            newFaces[newFaceI++] = face(newFace);
        }
        else
        {
            FatalErrorIn("shortEdgeFilter")
                << "Only " << newFace.size() << " in face " << faceI
                << exit(FatalError);
        }
    }

    newFaces.setSize(newFaceI);

    MeshedSurface<face> fMesh
    (
        xferMove(newPoints),
        xferMove(newFaces),
        xferCopy(List<surfZone>())
    );

    const Map<int>& fMeshPointMap = fMesh.meshPointMap();

    // Reset patchSizes_
    patchSizes_.clear();
    patchSizes_.setSize(patchNames_.size(), 0);

    label equalEdges = 0;
    label notFound = 0;
    label matches = 0;
    label negativeLabels = 0;

    forAll(newPointNumbers, pointI)
    {
        if (newPointNumbers[pointI] == -1)
        {
            WarningIn("shortEdgeFilter")
                << pointI << " will be deleted and " << newPointNumbers[pointI]
                << ", so it will not be replaced. "
                << "This will cause edges to be deleted." << endl;
        }
    }

    // Create new EdgeMap.
    Info<< "Creating new EdgeMap." << endl;
    EdgeMap<label> newMapEdgesRegion(mapEdgesRegion_.size());

    for
    (
        label bEdgeI = ms_.nInternalEdges();
        bEdgeI < edges.size();
        ++bEdgeI
    )
    {
        label p1 = meshPoints[edges[bEdgeI][0]];
        label p2 = meshPoints[edges[bEdgeI][1]];

        edge e(p1, p2);

        if (mapEdgesRegion_.found(e))
        {
            if
            (
                newPointNumbers[p1] != -1
             && newPointNumbers[p2] != -1
            )
            {
                if (newPointNumbers[p1] != newPointNumbers[p2])
                {
                    label region = mapEdgesRegion_.find(e)();
                    newMapEdgesRegion.insert
                    (
                        edge
                        (
                            fMeshPointMap[newPointNumbers[p1]],
                            fMeshPointMap[newPointNumbers[p2]]
                        ),
                        region
                    );
                    patchSizes_[region]++;
                    matches++;
                }
                else
                {
                    equalEdges++;
                }
            }
            else
            {
                negativeLabels++;
            }
        }
        else
        {
            notFound++;
        }
    }

    if (debug)
    {
        Info<< "EdgeMapping  :" << nl
            << "    Matches  : " << matches << nl
            << "    Equal    : " << equalEdges << nl
            << "    Negative : " << negativeLabels << nl
            << "    Not Found: " << notFound << endl;
    }

    mapEdgesRegion_.transfer(newMapEdgesRegion);

    ms_.transfer(fMesh);

    Info<< "    Maximum number of chained collapses = " << maxChain << endl;

    if (debug)
    {
        writeInfo(Info);
    }
}


void Foam::shortEdgeFilter2D::writeInfo(Ostream& os)
{
    os  << "Short Edge Filtering Information:" << nl
        << "    shortEdgeFilterFactor : " << shortEdgeFilterFactor_ << nl
        << "    edgeAttachedToBoundaryFactor : " << edgeAttachedToBoundaryFactor_
        << endl;

    forAll(patchNames_, patchI)
    {
        os  << "    Patch " << patchNames_[patchI]
            << ", size " << patchSizes_[patchI] << endl;
    }

    os  << "    There are " << mapEdgesRegion_.size()
        << " boundary edges." << endl;

    os  << "    Mesh Info:" << nl
        << "        Points:       " << ms_.nPoints() << nl
        << "        Faces:        " << ms_.size() << nl
        << "        Edges:        " << ms_.nEdges() << nl
        << "            Internal: " << ms_.nInternalEdges() << nl
        << "            External: " << ms_.nEdges() - ms_.nInternalEdges()
        << endl;
}


// ************************************************************************* //
