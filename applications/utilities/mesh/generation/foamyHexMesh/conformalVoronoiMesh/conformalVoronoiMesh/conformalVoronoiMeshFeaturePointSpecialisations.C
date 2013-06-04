/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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
) const
{
    List<extendedFeatureEdgeMesh::edgeStatus> allEdStat(pEds.size());

    forAll(pEds, i)
    {
        label edgeI = pEds[i];

        extendedFeatureEdgeMesh::edgeStatus& eS = allEdStat[i];

        eS = feMesh.getEdgeStatus(edgeI);

        pFEdgeTypes(eS)++;
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
    DynamicList<Vb>& pts
)
{
    if
    (
        !pFEdgesTypes.found(extendedFeatureEdgeMesh::EXTERNAL)
     || !pFEdgesTypes.found(extendedFeatureEdgeMesh::INTERNAL)
    )
    {
        return false;
    }

    if
    (
        pFEdgesTypes[extendedFeatureEdgeMesh::EXTERNAL] == 2
     && pFEdgesTypes[extendedFeatureEdgeMesh::INTERNAL] == 1
     && pEds.size() == 3
    )
    {
        Info<< "nExternal == 2 && nInternal == 1" << endl;

        const Foam::point& featPt = feMesh.points()[ptI];

        if (!positionOnThisProc(featPt))
        {
            return false;
        }

        label nVert = number_of_vertices();

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
        const Foam::point concaveEdgeLocalFeatPt =
            featPt + ppDist*concaveEdgeDir;

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

        pts.append
        (
            Vb(internalPtA, Vb::vtInternalFeaturePoint)
        );

        const Foam::point internalPtB =
            concaveEdgeExternalPt
          - 2.0*planeB.distance(concaveEdgeExternalPt)
            *concaveEdgePlaneBNormal;

        pts.append
        (
            Vb(internalPtB, Vb::vtInternalFeaturePoint)
        );

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
                    nVert--;
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

                        pts.append
                        (
                            Vb(externalPtD, Vb::vtExternalFeaturePoint)
                        );
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

                        pts.append
                        (
                            Vb(externalPtE, Vb::vtExternalFeaturePoint)
                        );
                    }
                }
            }
        }

        pts.append
        (
            Vb(concaveEdgeExternalPt, Vb::vtExternalFeaturePoint)
        );

        const scalar totalAngle = radToDeg
        (
            constant::mathematical::pi
          + radAngleBetween(concaveEdgePlaneANormal, concaveEdgePlaneBNormal)
        );

        if (totalAngle > foamyHexMeshControls().maxQuadAngle())
        {
            // Add additional mitreing points
            //scalar angleSign = 1.0;


            vector convexEdgesPlaneNormal =
                0.5*(convexEdgePlaneCNormal + convexEdgePlaneDNormal);

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

            pts.append
            (
                Vb(internalPtF, Vb::vtInternalFeaturePoint)
            );

            const Foam::point externalPtG =
                internalPtF
              + 2.0*planeM.distance(internalPtF)*convexEdgesPlaneNormal;

            pts.append
            (
                Vb(externalPtG, Vb::vtExternalFeaturePoint)
            );
        }

        if (debug)
        {
            for (label ptI = initialNumOfPoints; ptI < pts.size(); ++ptI)
            {
                Info<< "Point " << ptI << " : ";
                meshTools::writeOBJ(Info, topoint(pts[ptI].point()));
            }
        }

        return true;
    }
    else if
    (
        pFEdgesTypes[extendedFeatureEdgeMesh::EXTERNAL] == 1
     && pFEdgesTypes[extendedFeatureEdgeMesh::INTERNAL] == 2
     && pEds.size() == 3
    )
    {
        Info<< "nExternal == 1 && nInternal == 2" << endl;

        const Foam::point& featPt = feMesh.points()[ptI];

        if (!positionOnThisProc(featPt))
        {
            return false;
        }

        label nVert = number_of_vertices();

        const label initialNumOfPoints = pts.size();

        const scalar ppDist = pointPairDistance(featPt);

        const vectorField& normals = feMesh.normals();

        const labelListList& edgeNormals = feMesh.edgeNormals();

        label convexEdgeI = -1;
        labelList concaveEdgesI(2, -1);
        label nConcave = 0;

        forAll(pEds, i)
        {
            const extendedFeatureEdgeMesh::edgeStatus& eS = allEdStat[i];

            if (eS == extendedFeatureEdgeMesh::EXTERNAL)
            {
                convexEdgeI = pEds[i];
            }
            else if (eS == extendedFeatureEdgeMesh::INTERNAL)
            {
                concaveEdgesI[nConcave++] = pEds[i];
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

        const vector& convexEdgePlaneANormal =
            normals[edgeNormals[convexEdgeI][0]];

        const vector& convexEdgePlaneBNormal =
            normals[edgeNormals[convexEdgeI][1]];

        // Intersect planes parallel to the concave edge planes offset
        // by ppDist and the plane defined by featPt and the edge vector.
        plane planeA
        (
            featPt - ppDist*convexEdgePlaneANormal,
            convexEdgePlaneANormal
        );

        plane planeB
        (
            featPt - ppDist*convexEdgePlaneBNormal,
            convexEdgePlaneBNormal
        );

        const vector& convexEdgeDir = feMesh.edgeDirection
        (
            convexEdgeI,
            ptI
        );

        // Todo,needed later but want to get rid of this.
        const Foam::point convexEdgeLocalFeatPt =
            featPt + ppDist*convexEdgeDir;

        // Finding the nearest point on the intersecting line to the edge
        // point. Floating point errors often occur using planePlaneIntersect

        plane planeF(convexEdgeLocalFeatPt, convexEdgeDir);

        const Foam::point convexEdgeExternalPt = planeF.planePlaneIntersect
        (
            planeA,
            planeB
        );

        // Redefine planes to be on the feature surfaces to project through

        planeA = plane(featPt, convexEdgePlaneANormal);

        planeB = plane(featPt, convexEdgePlaneBNormal);

        const Foam::point internalPtA =
            convexEdgeExternalPt
          + 2.0*planeA.distance(convexEdgeExternalPt)
           *convexEdgePlaneANormal;

        pts.append
        (
            Vb(internalPtA, Vb::vtExternalFeaturePoint)
        );

        const Foam::point internalPtB =
            convexEdgeExternalPt
          + 2.0*planeB.distance(convexEdgeExternalPt)
           *convexEdgePlaneBNormal;

        pts.append
        (
            Vb(internalPtB, Vb::vtExternalFeaturePoint)
        );

        // Add the internal points

        Foam::point externalPtD;
        Foam::point externalPtE;

        vector concaveEdgePlaneCNormal(0,0,0);
        vector concaveEdgePlaneDNormal(0,0,0);

        const labelList& convexEdgeNormals = edgeNormals[convexEdgeI];
        const labelList& concaveEdgeANormals = edgeNormals[concaveEdgesI[0]];
        const labelList& concaveEdgeBNormals = edgeNormals[concaveEdgesI[1]];

        forAll(convexEdgeNormals, edgeNormalI)
        {
            bool concaveEdgeA = false;
            bool concaveEdgeB = false;

            forAll(concaveEdgeANormals, edgeAnormalI)
            {
                const vector& convexNormal
                    = normals[convexEdgeNormals[edgeNormalI]];
                const vector& concaveNormal
                    = normals[concaveEdgeANormals[edgeAnormalI]];

                Info<< "Angle between vectors = "
                    << degAngleBetween(convexNormal, concaveNormal) << endl;

                // Need a looser tolerance, because sometimes adjacent triangles
                // on the same surface will be slightly out of alignment.
                if (areParallel(convexNormal, concaveNormal, tolParallel))
                {
                    concaveEdgeA = true;
                }
            }

            forAll(concaveEdgeBNormals, edgeBnormalI)
            {
                const vector& convexNormal
                    = normals[convexEdgeNormals[edgeNormalI]];
                const vector& concaveNormal
                    = normals[concaveEdgeBNormals[edgeBnormalI]];

                Info<< "Angle between vectors = "
                    << degAngleBetween(convexNormal, concaveNormal) << endl;

                // Need a looser tolerance, because sometimes adjacent triangles
                // on the same surface will be slightly out of alignment.
                if (areParallel(convexNormal, concaveNormal, tolParallel))
                {
                    concaveEdgeB = true;
                }
            }

            if
            (
                (concaveEdgeA && concaveEdgeB)
             || (!concaveEdgeA && !concaveEdgeB)
            )
            {
                WarningIn
                (
                 "Foam::conformalVoronoiMesh"
                 "::createSpecialisedFeaturePoint"
                )   << "Both or neither of the concave edges share the convex "
                    << "edge's normal."
                    << " concaveEdgeA = " << concaveEdgeA
                    << " concaveEdgeB = " << concaveEdgeB
                    << endl;

                // Remove points that have just been added before returning
                for (label i = 0; i < 2; ++i)
                {
                    pts.remove();
                    nVert--;
                }

                return false;
            }

            if (concaveEdgeA)
            {
                forAll(concaveEdgeANormals, edgeAnormalI)
                {
                    const vector& convexNormal
                        = normals[convexEdgeNormals[edgeNormalI]];
                    const vector& concaveNormal
                        = normals[concaveEdgeANormals[edgeAnormalI]];

                    if
                    (
                        !areParallel(convexNormal, concaveNormal, tolParallel)
                    )
                    {
                        concaveEdgePlaneCNormal = concaveNormal;

                        plane planeC(featPt, concaveEdgePlaneCNormal);

                        externalPtD =
                            internalPtA
                          - 2.0*planeC.distance(internalPtA)
                           *concaveEdgePlaneCNormal;

                        pts.append
                        (
                            Vb(externalPtD, Vb::vtInternalFeaturePoint)
                        );
                    }
                }
            }

            if (concaveEdgeB)
            {
                forAll(concaveEdgeBNormals, edgeBnormalI)
                {
                    const vector& convexNormal
                        = normals[convexEdgeNormals[edgeNormalI]];
                    const vector& concaveNormal
                        = normals[concaveEdgeBNormals[edgeBnormalI]];

                    if
                    (
                        !areParallel(convexNormal, concaveNormal, tolParallel)
                    )
                    {
                        concaveEdgePlaneDNormal = concaveNormal;

                        plane planeD(featPt, concaveEdgePlaneDNormal);

                        externalPtE =
                            internalPtB
                          - 2.0*planeD.distance(internalPtB)
                           *concaveEdgePlaneDNormal;

                        pts.append
                        (
                            Vb(externalPtE, Vb::vtInternalFeaturePoint)
                        );
                    }
                }
            }
        }

        pts.append
        (
            Vb(convexEdgeExternalPt, Vb::vtInternalFeaturePoint)
        );

        const scalar totalAngle = radToDeg
        (
            constant::mathematical::pi
          + radAngleBetween(convexEdgePlaneANormal, convexEdgePlaneBNormal)
        );

        if (totalAngle > foamyHexMeshControls().maxQuadAngle())
        {
            // Add additional mitreing points
            //scalar angleSign = 1.0;


            vector convexEdgesPlaneNormal =
                0.5*(concaveEdgePlaneCNormal + concaveEdgePlaneDNormal);

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
                convexEdgeExternalPt
            //+ (2.0 + guard)*(concaveEdgeLocalFeatPt - concaveEdgeExternalPt);
              + 2.0*(convexEdgeLocalFeatPt - convexEdgeExternalPt);

            pts.append
            (
                Vb(internalPtF, Vb::vtExternalFeaturePoint)
            );

            const Foam::point externalPtG =
                internalPtF
              - 2.0*planeM.distance(internalPtF)*convexEdgesPlaneNormal;

            pts.append
            (
                Vb(externalPtG, Vb::vtInternalFeaturePoint)
            );
        }

        if (debug)
        {
            for (label ptI = initialNumOfPoints; ptI < pts.size(); ++ptI)
            {
                Info<< "Point " << ptI << " "
                    << indexedVertexEnum::vertexTypeNames_[pts[ptI].type()]
                    << " : ";
                meshTools::writeOBJ(Info, topoint(pts[ptI].point()));
            }
        }

        return true;
    }

    return false;
}


// ************************************************************************* //
