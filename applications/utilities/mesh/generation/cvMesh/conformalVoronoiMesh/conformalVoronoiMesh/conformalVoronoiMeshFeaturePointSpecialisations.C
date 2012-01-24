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

\*---------------------------------------------------------------------------*/

#include "conformalVoronoiMesh.H"
#include "vectorTools.H"

using namespace Foam::vectorTools;

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::List<Foam::extendedFeatureEdgeMesh::edgeStatus>
Foam::conformalVoronoiMesh::calcPointFeatureEdgesTypes
(
    const extendedFeatureEdgeMesh& feMesh,
    const labelList& pEds,
    pointFeatureEdgesTypes& pFEdgeTypes
)
{
    List<extendedFeatureEdgeMesh::edgeStatus> allEdStat(pEds.size());

    forAll(pEds, i)
    {
        label edgeI = pEds[i];

        extendedFeatureEdgeMesh::edgeStatus& eS = allEdStat[i];

        eS = feMesh.getEdgeStatus(edgeI);

        switch (eS)
        {
            case extendedFeatureEdgeMesh::EXTERNAL:
            {
                pFEdgeTypes.nExternal++;
                break;
            }
            case extendedFeatureEdgeMesh::INTERNAL:
            {
                pFEdgeTypes.nInternal++;
                break;
            }
            case extendedFeatureEdgeMesh::FLAT:
            {
                pFEdgeTypes.nFlat++;
                break;
            }
            case extendedFeatureEdgeMesh::OPEN:
            {
                pFEdgeTypes.nOpen++;
                break;
            }
            case extendedFeatureEdgeMesh::MULTIPLE:
            {
                pFEdgeTypes.nMultiple++;
                break;
            }
            case extendedFeatureEdgeMesh::NONE:
            {
                pFEdgeTypes.nNonFeature++;
                break;
            }
        }
    }

    if (debug)
    {
        Info<< pFEdgeTypes << endl;
    }

    return allEdStat;
}


bool Foam::conformalVoronoiMesh::createSpecialisedFeaturePoint
(
    const extendedFeatureEdgeMesh& feMesh,
    const labelList& pEds,
    const pointFeatureEdgesTypes& pFEdgesTypes,
    const List<extendedFeatureEdgeMesh::edgeStatus>& allEdStat,
    const label ptI,
    DynamicList<Foam::point>& pts,
    DynamicList<label>& indices,
    DynamicList<label>& types
)
{
    if
    (
        pFEdgesTypes.nExternal == 2
     && pFEdgesTypes.nInternal == 1
     && pEds.size() == 3
    )
    {
        Info<< "nExternal == 2 && nInternal == 1" << endl;

        const Foam::point& featPt = feMesh.points()[ptI];

        if (!positionOnThisProc(featPt))
        {
            return false;
        }

        const label initialNumOfPoints = pts.size();

        const scalar ppDist = pointPairDistance(featPt);

        const vectorField& normals = feMesh.normals();

        const labelListList& edgeNormals = feMesh.edgeNormals();

        label concaveEdgeI = -1;
        labelList convexEdgesI(2, -1);
        label nConvex = 0;

        forAll(pEds, i)
        {
            const extendedFeatureEdgeMesh::edgeStatus& eS = allEdStat[i];

            if (eS == extendedFeatureEdgeMesh::INTERNAL)
            {
                concaveEdgeI = pEds[i];
            }
            else if (eS == extendedFeatureEdgeMesh::EXTERNAL)
            {
                convexEdgesI[nConvex++] = pEds[i];
            }
            else if (eS == extendedFeatureEdgeMesh::FLAT)
            {
                WarningIn("Foam::conformalVoronoiMesh::"
                    "createSpecialisedFeaturePoint")
                    << "Edge " << eS << " is flat"
                    << endl;
            }
            else
            {
                FatalErrorIn("Foam::conformalVoronoiMesh::"
                    "createSpecialisedFeaturePoint")
                    << "Edge " << eS << " not concave/convex"
                    << exit(FatalError);
            }
        }

        const vector& concaveEdgePlaneANormal =
            normals[edgeNormals[concaveEdgeI][0]];

        const vector& concaveEdgePlaneBNormal =
            normals[edgeNormals[concaveEdgeI][1]];

        // Intersect planes parallel to the concave edge planes offset
        // by ppDist and the plane defined by featPt and the edge vector.
        plane planeA
        (
            featPt + ppDist*concaveEdgePlaneANormal,
            concaveEdgePlaneANormal
        );

        plane planeB
        (
            featPt + ppDist*concaveEdgePlaneBNormal,
            concaveEdgePlaneBNormal
        );

        const vector& concaveEdgeDir = feMesh.edgeDirection
        (
            concaveEdgeI,
            ptI
        );

        // Todo,needed later but want to get rid of this.
        const Foam::point concaveEdgeLocalFeatPt = featPt + ppDist*concaveEdgeDir;

        // Finding the nearest point on the intersecting line to the edge
        // point. Floating point errors often occur using planePlaneIntersect

        plane planeF(concaveEdgeLocalFeatPt, concaveEdgeDir);

        const Foam::point concaveEdgeExternalPt = planeF.planePlaneIntersect
        (
            planeA,
            planeB
        );

        // Redefine planes to be on the feature surfaces to project through

        planeA = plane(featPt, concaveEdgePlaneANormal);

        planeB = plane(featPt, concaveEdgePlaneBNormal);

        const Foam::point internalPtA =
            concaveEdgeExternalPt
          - 2.0*planeA.distance(concaveEdgeExternalPt)
            *concaveEdgePlaneANormal;

        pts.append(internalPtA);
        indices.append(0);
        types.append(4);

        const Foam::point internalPtB =
            concaveEdgeExternalPt
          - 2.0*planeB.distance(concaveEdgeExternalPt)
            *concaveEdgePlaneBNormal;

        pts.append(internalPtB);
        indices.append(0);
        types.append(3);

        // Add the external points

        Foam::point externalPtD;
        Foam::point externalPtE;

        vector convexEdgePlaneCNormal(0,0,0);
        vector convexEdgePlaneDNormal(0,0,0);

        const labelList& concaveEdgeNormals = edgeNormals[concaveEdgeI];
        const labelList& convexEdgeANormals = edgeNormals[convexEdgesI[0]];
        const labelList& convexEdgeBNormals = edgeNormals[convexEdgesI[1]];

        forAll(concaveEdgeNormals, edgeNormalI)
        {
            bool convexEdgeA = false;
            bool convexEdgeB = false;

            forAll(convexEdgeANormals, edgeAnormalI)
            {
                const vector& concaveNormal
                    = normals[concaveEdgeNormals[edgeNormalI]];
                const vector& convexNormal
                    = normals[convexEdgeANormals[edgeAnormalI]];

                Info<< "Angle between vectors = "
                    << degAngleBetween(concaveNormal, convexNormal) << endl;

                // Need a looser tolerance, because sometimes adjacent triangles
                // on the same surface will be slightly out of alignment.
                if (areParallel(concaveNormal, convexNormal, tolParallel))
                {
                    convexEdgeA = true;
                }
            }

            forAll(convexEdgeBNormals, edgeBnormalI)
            {
                const vector& concaveNormal
                    = normals[concaveEdgeNormals[edgeNormalI]];
                const vector& convexNormal
                    = normals[convexEdgeBNormals[edgeBnormalI]];

                Info<< "Angle between vectors = "
                    << degAngleBetween(concaveNormal, convexNormal) << endl;

                // Need a looser tolerance, because sometimes adjacent triangles
                // on the same surface will be slightly out of alignment.
                if (areParallel(concaveNormal, convexNormal, tolParallel))
                {
                    convexEdgeB = true;
                }
            }

            if ((convexEdgeA && convexEdgeB) || (!convexEdgeA && !convexEdgeB))
            {
                WarningIn
                    (
                     "Foam::conformalVoronoiMesh"
                     "::createSpecialisedFeaturePoint"
                    )
                    << "Both or neither of the convex edges share the concave "
                    << "edge's normal."
                    << " convexEdgeA = " << convexEdgeA
                    << " convexEdgeB = " << convexEdgeB
                    << endl;

                // Remove points that have just been added before returning
                for (label i = 0; i < 2; ++i)
                {
                    pts.remove();
                    indices.remove();
                    types.remove();
                }

                return false;
            }

            if (convexEdgeA)
            {
                forAll(convexEdgeANormals, edgeAnormalI)
                {
                    const vector& concaveNormal
                        = normals[concaveEdgeNormals[edgeNormalI]];
                    const vector& convexNormal
                        = normals[convexEdgeANormals[edgeAnormalI]];

                    if
                    (
                        !areParallel(concaveNormal, convexNormal, tolParallel)
                    )
                    {
                        convexEdgePlaneCNormal = convexNormal;

                        plane planeC(featPt, convexEdgePlaneCNormal);

                        externalPtD =
                            internalPtA
                          + 2.0*planeC.distance(internalPtA)
                           *convexEdgePlaneCNormal;

                        pts.append(externalPtD);
                        indices.append(0);
                        types.append(-2);
                    }
                }
            }

            if (convexEdgeB)
            {
                forAll(convexEdgeBNormals, edgeBnormalI)
                {
                    const vector& concaveNormal
                        = normals[concaveEdgeNormals[edgeNormalI]];
                    const vector& convexNormal
                        = normals[convexEdgeBNormals[edgeBnormalI]];

                    if
                    (
                        !areParallel(concaveNormal, convexNormal, tolParallel)
                    )
                    {
                        convexEdgePlaneDNormal = convexNormal;

                        plane planeD(featPt, convexEdgePlaneDNormal);

                        externalPtE =
                            internalPtB
                          + 2.0*planeD.distance(internalPtB)
                           *convexEdgePlaneDNormal;

                        pts.append(externalPtE);
                        indices.append(0);
                        types.append(-2);
                    }
                }
            }
        }

        pts.append(concaveEdgeExternalPt);
        indices.append(0);
        types.append(-4);

        const scalar totalAngle = radToDeg
        (
            constant::mathematical::pi
          + radAngleBetween(concaveEdgePlaneANormal, concaveEdgePlaneBNormal)
        );

        if (totalAngle > cvMeshControls().maxQuadAngle())
        {
            // Add additional mitering points
            //scalar angleSign = 1.0;

            Info<<convexEdgePlaneCNormal << " " <<convexEdgePlaneDNormal << endl;

            vector convexEdgesPlaneNormal = 0.5*(convexEdgePlaneCNormal + convexEdgePlaneDNormal);
            plane planeM(featPt, convexEdgesPlaneNormal);

//            if
//            (
//                geometryToConformTo_.outside
//                (
//                    featPt - convexEdgesPlaneNormal*ppDist
//                )
//            )
//            {
//                angleSign = -1.0;
//            }

//            scalar phi =
//                angleSign*acos(concaveEdgeDir & -convexEdgesPlaneNormal);
//
//            scalar guard =
//            (
//                1.0 + sin(phi)*ppDist/mag
//                (
//                    concaveEdgeLocalFeatPt - concaveEdgeExternalPt
//                )
//            )/cos(phi) - 1.0;

            const Foam::point internalPtF =
                concaveEdgeExternalPt
            //+ (2.0 + guard)*(concaveEdgeLocalFeatPt - concaveEdgeExternalPt);
              + 2.0*(concaveEdgeLocalFeatPt - concaveEdgeExternalPt);

            pts.append(internalPtF);
            indices.append(0);
            types.append(1);

            const Foam::point externalPtG =
                internalPtF
              + 2.0*planeM.distance(internalPtF)*convexEdgesPlaneNormal;

            pts.append(externalPtG);
            indices.append(0);
            types.append(-1);
        }

        if (debug)
        {
            for (label ptI = initialNumOfPoints; ptI < pts.size(); ++ptI)
            {
                Info<< "Point " << ptI << " : ";
                meshTools::writeOBJ(Info, pts[ptI]);
            }
        }

        return true;
    }
    else if (pFEdgesTypes.nExternal == 1 && pFEdgesTypes.nInternal == 2)
    {
        // Info<< "nExternal == 1 && nInternal == 2" << endl;

        return false;
    }

    return false;
}


// ************************************************************************* //
