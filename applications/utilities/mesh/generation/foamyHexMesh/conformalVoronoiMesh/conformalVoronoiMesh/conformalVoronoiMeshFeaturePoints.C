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
#include "triangle.H"
#include "tetrahedron.H"
#include "const_circulator.H"

using namespace Foam::vectorTools;

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::conformalVoronoiMesh::createEdgePointGroup
(
    const extendedFeatureEdgeMesh& feMesh,
    const pointIndexHit& edHit,
    DynamicList<Vb>& pts
)
{
    if (foamyHexMeshControls().circulateEdges())
    {
        createEdgePointGroupByCirculating(feMesh, edHit, pts);
    }
    else
    {
        label edgeI = edHit.index();

        extendedFeatureEdgeMesh::edgeStatus edStatus =
            feMesh.getEdgeStatus(edgeI);

        switch (edStatus)
        {
            case extendedFeatureEdgeMesh::EXTERNAL:
            {
                createExternalEdgePointGroup(feMesh, edHit, pts);
                break;
            }
            case extendedFeatureEdgeMesh::INTERNAL:
            {
                createInternalEdgePointGroup(feMesh, edHit, pts);
                break;
            }
            case extendedFeatureEdgeMesh::FLAT:
            {
                createFlatEdgePointGroup(feMesh, edHit, pts);
                break;
            }
            case extendedFeatureEdgeMesh::OPEN:
            {
                createOpenEdgePointGroup(feMesh, edHit, pts);
                break;
            }
            case extendedFeatureEdgeMesh::MULTIPLE:
            {
                createMultipleEdgePointGroup(feMesh, edHit, pts);
                break;
            }
            case extendedFeatureEdgeMesh::NONE:
            {
                break;
            }
        }
    }
}


bool Foam::conformalVoronoiMesh::meshableRegion
(
    const plane::side side,
    const extendedFeatureEdgeMesh::sideVolumeType volType
) const
{
    switch (volType)
    {
        case extendedFeatureEdgeMesh::INSIDE:
        {
            return (side == plane::FLIP) ? true : false;
        }
        case extendedFeatureEdgeMesh::OUTSIDE:
        {
            return (side == plane::NORMAL) ? true : false;
        }
        case extendedFeatureEdgeMesh::BOTH:
        {
            return true;
        }
        case extendedFeatureEdgeMesh::NEITHER:
        {
            return false;
        }
    }
}


bool Foam::conformalVoronoiMesh::regionIsInside
(
    const extendedFeatureEdgeMesh::sideVolumeType volTypeA,
    const vector& normalA,
    const extendedFeatureEdgeMesh::sideVolumeType volTypeB,
    const vector& normalB,
    const vector& masterPtVec
) const
{
    plane::side sideA
    (
        ((masterPtVec & normalA) <= 0) ? plane::FLIP : plane::NORMAL
    );

    plane::side sideB
    (
        ((masterPtVec & normalB) <= 0) ? plane::FLIP : plane::NORMAL
    );

    const bool meshableRegionA = meshableRegion(sideA, volTypeA);
    const bool meshableRegionB = meshableRegion(sideB, volTypeB);

    if (meshableRegionA == meshableRegionB)
    {
        return meshableRegionA;
    }
    else
    {
        WarningIn
        (
            "Foam::conformalVoronoiMesh::regionIsInside"
            "(volTypeA, normalA, volTypeB, normalB, masterPtVec)"
        )   << ""
            << endl;

        return false;
    }
}


void Foam::conformalVoronoiMesh::createEdgePointGroupByCirculating
(
    const extendedFeatureEdgeMesh& feMesh,
    const pointIndexHit& edHit,
    DynamicList<Vb>& pts
)
{
    typedef Foam::indexedVertexEnum::vertexType         vertexType;
    typedef extendedFeatureEdgeMesh::sideVolumeType     sideVolumeType;

    const Foam::point& edgePt = edHit.hitPoint();
    const label edgeI = edHit.index();

    scalar ppDist = pointPairDistance(edgePt);

    const vectorField& feNormals = feMesh.normals();
    const vector& edDir = feMesh.edgeDirections()[edgeI];
    const labelList& edNormalIs = feMesh.edgeNormals()[edgeI];
    const labelList& feNormalDirections = feMesh.normalDirections()[edgeI];

    const PackedList<2>& normalVolumeTypes = feMesh.normalVolumeTypes();

    const_circulator<labelList> circ(edNormalIs);
    const_circulator<labelList> circNormalDirs(feNormalDirections);

    Map<Foam::point> masterPoints;
    Map<vertexType> masterPointsTypes;
    Map<plane> masterPointReflectionsPrev;
    Map<plane> masterPointReflectionsNext;

//    Info<< "Edge = " << edHit << ", edDir = " << edDir << endl;
//    Info<< "    edNorms = " << edNormalIs.size() << endl;

    bool addedMasterPreviously = false;
    label initialRegion = -1;

    if (circ.size()) do
    {
        const sideVolumeType volType =
            sideVolumeType(normalVolumeTypes[circ()]);
        const sideVolumeType nextVolType =
            sideVolumeType(normalVolumeTypes[circ.next()]);

        const vector& normal = feNormals[circ()];
        const vector& nextNormal = feNormals[circ.next()];

        vector normalDir = (normal ^ edDir);
        normalDir *= circNormalDirs()/mag(normalDir);

        vector nextNormalDir = (nextNormal ^ edDir);
        nextNormalDir *= circNormalDirs.next()/mag(nextNormalDir);

//        Info<< "    " << circ() << " " << circ.next() << nl
//            << "    " << circNormalDirs() << " " << circNormalDirs.next()
//            << nl
//            << "    normals = " << normal << " " << nextNormal << nl
//            << " normalDirs = " << normalDir << " " << nextNormalDir << nl
//            << "      cross = " << (normalDir ^ nextNormalDir) << nl
//            << "    "
//            << extendedFeatureEdgeMesh::sideVolumeTypeNames_[volType] << " "
//            << extendedFeatureEdgeMesh::sideVolumeTypeNames_[nextVolType]
//            << endl;

        // Calculate master point
        vector masterPtVec(normalDir + nextNormalDir);
        masterPtVec /= mag(masterPtVec) + SMALL;

        if (((normalDir ^ nextNormalDir) & edDir) < SMALL)
        {
//            Info<< "    IGNORE REGION" << endl;
            addedMasterPreviously = false;

            if (circ.size() == 2 && mag((normal & nextNormal) - 1) < SMALL)
            {
                // Add an extra point

                const vector n = 0.5*(normal + nextNormal);

                vector s = ppDist*(edDir ^ n);

                createPointPair(ppDist, edgePt + s, n, pts);
                createPointPair(ppDist, edgePt - s, n, pts);

                break;
            }

            continue;
        }

        if (mag(masterPtVec) < SMALL)
        {
            if (circ.size() == 2)
            {
                // Add an extra point

                // Average normal to remove any bias to one side, although as
                // it is a flat edge, the normals should be essentially the same
                const vector n = 0.5*(normal + nextNormal);

                // Direction along the surface to the control point, sense of
                // direction not important, as +s and -s can be used because it
                // is a flat edge
                vector s = ppDist*(edDir ^ n);

                createPointPair(ppDist, edgePt + s, n, pts);
                createPointPair(ppDist, edgePt - s, n, pts);

                break;
            }
            else
            {
//                Info<< "    IGNORE REGION" << endl;
                addedMasterPreviously = false;
                continue;
            }
        }

        const Foam::point masterPt = edgePt + ppDist*masterPtVec;

        // Check that region is inside or outside
        const bool inside =
            regionIsInside
            (
                volType,
                normal,
                nextVolType,
                nextNormal,
                masterPtVec
            );

        // Specialise for size = 1 && baffle
        if (mag((normalDir & nextNormalDir) - 1) < SMALL)
        {
            if (inside)
            {
//                Info<< "Specialise for size 1 and baffle" << endl;

                vector s = ppDist*(edDir ^ normal);

                plane facePlane(edgePt, normal);

                Foam::point pt1 = edgePt + s + ppDist*normal;
                Foam::point pt2 = edgePt - s + ppDist*normal;

                Foam::point pt3 = reflectPointInPlane(pt1, facePlane);
                Foam::point pt4 = reflectPointInPlane(pt2, facePlane);

                pts.append(Vb(pt1, Vb::vtInternalFeatureEdge));
                pts.append(Vb(pt2, Vb::vtInternalFeatureEdge));
                pts.append(Vb(pt3, Vb::vtInternalFeatureEdge));
                pts.append(Vb(pt4, Vb::vtInternalFeatureEdge));

                break;
            }
            else
            {
                WarningIn
                (
                    "Foam::conformalVoronoiMesh::"
                    "createEdgePointGroupByCirculating"
                    "("
                    "    const extendedFeatureEdgeMesh&,"
                    "    const pointIndexHit&,"
                    "    DynamicList<Vb>&"
                    ")"
                )   << "Faces are parallel but master point is not inside"
                    << endl;
            }
        }

        if (!addedMasterPreviously)
        {
            if (initialRegion == -1)
            {
                initialRegion = circ.nRotations();
            }

            addedMasterPreviously = true;

            masterPoints.insert(circ.nRotations(), masterPt);
            masterPointsTypes.insert
            (
                circ.nRotations(),
                inside
              ? Vb::vtInternalFeatureEdge
              : Vb::vtExternalFeatureEdge
            );

            masterPointReflectionsPrev.insert
            (
                circ.nRotations(),
                plane(edgePt, normal)
            );

            masterPointReflectionsNext.insert
            (
                circ.nRotations(),
                plane(edgePt, nextNormal)
            );
        }
        else if (addedMasterPreviously)
        {
            addedMasterPreviously = true;

            masterPointReflectionsNext.erase(circ.nRotations() - 1);

            // Shift the master point to be normal to the plane between it and
            // the previous master point
            // Should be the intersection of the normal and the plane with the
            // new master point in it.

            plane p(masterPoints[circ.nRotations() - 1], normalDir);
            plane::ray r(edgePt, masterPt - edgePt);

            scalar cutPoint = p.normalIntersect(r);

            masterPoints.insert
            (
                circ.nRotations(),
                edgePt + cutPoint*(masterPt - edgePt)
            );

            masterPointsTypes.insert
            (
                circ.nRotations(),
                inside
              ? Vb::vtInternalFeatureEdge
              : Vb::vtExternalFeatureEdge
            );

            masterPointReflectionsNext.insert
            (
                circ.nRotations(),
                plane(edgePt, nextNormal)
            );
        }

        if
        (
            masterPoints.size() > 1
         && inside
         && circ.nRotations() == circ.size() - 1
        )
        {
            if (initialRegion == 0)
            {
                plane p(masterPoints[initialRegion], nextNormalDir);
                plane::ray r(edgePt, masterPt - edgePt);

                scalar cutPoint = p.normalIntersect(r);

                masterPoints[circ.nRotations()] =
                    edgePt + cutPoint*(masterPt - edgePt);

                // Remove the first reflection plane if we are no longer
                // circulating
                masterPointReflectionsPrev.erase(initialRegion);
                masterPointReflectionsNext.erase(circ.nRotations());
            }
            else
            {

            }
        }
    }
    while
    (
        circ.circulate(CirculatorBase::CLOCKWISE),
        circNormalDirs.circulate(CirculatorBase::CLOCKWISE)
    );


    forAllConstIter(Map<Foam::point>, masterPoints, iter)
    {
        const Foam::point& pt = masterPoints[iter.key()];
        const vertexType ptType = masterPointsTypes[iter.key()];

//        Info<< "    Adding Master " << iter.key() << " " << pt << " "
//            << indexedVertexEnum::vertexTypeNames_[ptType] << endl;

        pts.append(Vb(pt, ptType));

        const vertexType reflectedPtType =
        (
            ptType == Vb::vtInternalFeatureEdge
          ? Vb::vtExternalFeatureEdge
          : Vb::vtInternalFeatureEdge
        );

        if (masterPointReflectionsPrev.found(iter.key()))
        {
            const Foam::point reflectedPt =
                reflectPointInPlane(pt, masterPointReflectionsPrev[iter.key()]);

//            Info<< "        Adding Prev " << reflectedPt << " "
//                << indexedVertexEnum::vertexTypeNames_[reflectedPtType]
//                << endl;

            pts.append(Vb(reflectedPt, reflectedPtType));
        }

        if (masterPointReflectionsNext.found(iter.key()))
        {
            const Foam::point reflectedPt =
                reflectPointInPlane(pt, masterPointReflectionsNext[iter.key()]);

//            Info<< "        Adding Next " << reflectedPt << " "
//                << indexedVertexEnum::vertexTypeNames_[reflectedPtType]
//                << endl;

            pts.append(Vb(reflectedPt, reflectedPtType));
        }
    }

//    pts.append(Vb(edgePt, Vb::vtExternalFeatureEdge));
}


void Foam::conformalVoronoiMesh::createExternalEdgePointGroup
(
    const extendedFeatureEdgeMesh& feMesh,
    const pointIndexHit& edHit,
    DynamicList<Vb>& pts
)
{
    const Foam::point& edgePt = edHit.hitPoint();

    scalar ppDist = pointPairDistance(edgePt);

    const vectorField& feNormals = feMesh.normals();
    const labelList& edNormalIs = feMesh.edgeNormals()[edHit.index()];

    // As this is an external edge, there are two normals by definition
    const vector& nA = feNormals[edNormalIs[0]];
    const vector& nB = feNormals[edNormalIs[1]];

    if (areParallel(nA, nB))
    {
        // The normals are nearly parallel, so this is too sharp a feature to
        // conform to.
        return;
    }

    // Normalised distance of reference point from edge point
    vector refVec((nA + nB)/(1 + (nA & nB)));

    if (magSqr(refVec) > sqr(5.0))
    {
        // Limit the size of the conformation
        ppDist *= 5.0/mag(refVec);

        // Pout<< nl << "createExternalEdgePointGroup limit "
        //     << "edgePt " << edgePt << nl
        //     << "refVec " << refVec << nl
        //     << "mag(refVec) " << mag(refVec) << nl
        //     << "ppDist " << ppDist << nl
        //     << "nA " << nA << nl
        //     << "nB " << nB << nl
        //     << "(nA & nB) " << (nA & nB) << nl
        //     << endl;
    }

    // Convex. So refPt will be inside domain and hence a master point
    Foam::point refPt = edgePt - ppDist*refVec;

    // Result when the points are eventually inserted.
    // Add number_of_vertices() at insertion of first vertex to all numbers:
    // pt           index type
    // refPt        0     1
    // reflectedA   1     0
    // reflectedB   2     0

    // Insert the master point pairing the the first slave

    if (!geometryToConformTo_.inside(refPt))
    {
        return;
    }

    pts.append
    (
        Vb
        (
            refPt,
            vertexCount() + pts.size(),
            Vb::vtInternalFeatureEdge,
            Pstream::myProcNo()
        )
    );

    // Insert the slave points by reflecting refPt in both faces.
    // with each slave refering to the master

    Foam::point reflectedA = refPt + 2*ppDist*nA;
    pts.append
    (
        Vb
        (
            reflectedA,
            vertexCount() + pts.size(),
            Vb::vtExternalFeatureEdge,
            Pstream::myProcNo()
        )
    );

    Foam::point reflectedB = refPt + 2*ppDist*nB;
    pts.append
    (
        Vb
        (
            reflectedB,
            vertexCount() + pts.size(),
            Vb::vtExternalFeatureEdge,
            Pstream::myProcNo()
        )
    );
}


void Foam::conformalVoronoiMesh::createInternalEdgePointGroup
(
    const extendedFeatureEdgeMesh& feMesh,
    const pointIndexHit& edHit,
    DynamicList<Vb>& pts
)
{
    const Foam::point& edgePt = edHit.hitPoint();

    scalar ppDist = pointPairDistance(edgePt);

    const vectorField& feNormals = feMesh.normals();
    const labelList& edNormalIs = feMesh.edgeNormals()[edHit.index()];

    // As this is an external edge, there are two normals by definition
    const vector& nA = feNormals[edNormalIs[0]];
    const vector& nB = feNormals[edNormalIs[1]];

    if (areParallel(nA, nB))
    {
        // The normals are nearly parallel, so this is too sharp a feature to
        // conform to.
        return;
    }

    // Normalised distance of reference point from edge point
    vector refVec((nA + nB)/(1 + (nA & nB)));

    if (magSqr(refVec) > sqr(5.0))
    {
        // Limit the size of the conformation
        ppDist *= 5.0/mag(refVec);

        // Pout<< nl << "createInternalEdgePointGroup limit "
        //     << "edgePt " << edgePt << nl
        //     << "refVec " << refVec << nl
        //     << "mag(refVec) " << mag(refVec) << nl
        //     << "ppDist " << ppDist << nl
        //     << "nA " << nA << nl
        //     << "nB " << nB << nl
        //     << "(nA & nB) " << (nA & nB) << nl
        //     << endl;
    }

    // Concave. master and reflected points inside the domain.
    Foam::point refPt = edgePt - ppDist*refVec;

    // Generate reflected master to be outside.
    Foam::point reflMasterPt = refPt + 2*(edgePt - refPt);

    // Reflect reflMasterPt in both faces.
    Foam::point reflectedA = reflMasterPt - 2*ppDist*nA;

    Foam::point reflectedB = reflMasterPt - 2*ppDist*nB;

    scalar totalAngle =
        radToDeg(constant::mathematical::pi + radAngleBetween(nA, nB));

    // Number of quadrants the angle should be split into
    int nQuads = int(totalAngle/foamyHexMeshControls().maxQuadAngle()) + 1;

    // The number of additional master points needed to obtain the
    // required number of quadrants.
    int nAddPoints = min(max(nQuads - 2, 0), 2);

    // Add number_of_vertices() at insertion of first vertex to all numbers:
    // Result for nAddPoints 1 when the points are eventually inserted
    // pt           index type
    // reflectedA   0     2
    // reflectedB   1     2
    // reflMasterPt 2     0

    // Result for nAddPoints 1 when the points are eventually inserted
    // pt           index type
    // reflectedA   0     3
    // reflectedB   1     3
    // refPt        2     3
    // reflMasterPt 3     0

    // Result for nAddPoints 2 when the points are eventually inserted
    // pt           index type
    // reflectedA   0     4
    // reflectedB   1     4
    // reflectedAa  2     4
    // reflectedBb  3     4
    // reflMasterPt 4     0

    if
    (
        !geometryToConformTo_.inside(reflectedA)
     || !geometryToConformTo_.inside(reflectedB)
    )
    {
        return;
    }

    // Master A is inside.
    pts.append
    (
        Vb
        (
            reflectedA,
            vertexCount() + pts.size(),
            Vb::vtInternalFeatureEdge,
            Pstream::myProcNo()
        )
    );

    // Master B is inside.
    pts.append
    (
        Vb
        (
            reflectedB,
            vertexCount() + pts.size(),
            Vb::vtInternalFeatureEdge,
            Pstream::myProcNo()
        )
    );

    // Slave is outside.
    pts.append
    (
        Vb
        (
            reflMasterPt,
            vertexCount() + pts.size(),
            Vb::vtExternalFeatureEdge,
            Pstream::myProcNo()
        )
    );

    if (nAddPoints == 1)
    {
        // One additinal point is the reflection of the slave point,
        // i.e. the original reference point
        pts.append
        (
            Vb
            (
                refPt,
                vertexCount() + pts.size(),
                Vb::vtInternalFeatureEdge,
                Pstream::myProcNo()
            )
        );
    }
    else if (nAddPoints == 2)
    {
        Foam::point reflectedAa = refPt + ppDist*nB;
        pts.append
        (
            Vb
            (
                reflectedAa,
                vertexCount() + pts.size(),
                Vb::vtInternalFeatureEdge,
                Pstream::myProcNo()
            )
        );

        Foam::point reflectedBb = refPt + ppDist*nA;
        pts.append
        (
            Vb
            (
                reflectedBb,
                vertexCount() + pts.size(),
                Vb::vtInternalFeatureEdge,
                Pstream::myProcNo()
            )
        );
    }
}


void Foam::conformalVoronoiMesh::createFlatEdgePointGroup
(
    const extendedFeatureEdgeMesh& feMesh,
    const pointIndexHit& edHit,
    DynamicList<Vb>& pts
)
{
    const Foam::point& edgePt = edHit.hitPoint();

    const scalar ppDist = pointPairDistance(edgePt);

    const vectorField& feNormals = feMesh.normals();
    const labelList& edNormalIs = feMesh.edgeNormals()[edHit.index()];

    // As this is a flat edge, there are two normals by definition
    const vector& nA = feNormals[edNormalIs[0]];
    const vector& nB = feNormals[edNormalIs[1]];

    // Average normal to remove any bias to one side, although as this
    // is a flat edge, the normals should be essentially the same
    const vector n = 0.5*(nA + nB);

    // Direction along the surface to the control point, sense of edge
    // direction not important, as +s and -s can be used because this
    // is a flat edge
    vector s = ppDist*(feMesh.edgeDirections()[edHit.index()] ^ n);

    createPointPair(ppDist, edgePt + s, n, pts);
    createPointPair(ppDist, edgePt - s, n, pts);
}


void Foam::conformalVoronoiMesh::createOpenEdgePointGroup
(
    const extendedFeatureEdgeMesh& feMesh,
    const pointIndexHit& edHit,
    DynamicList<Vb>& pts
)
{
    // Assume it is a baffle and insert flat edge point pairs
    const Foam::point& edgePt = edHit.hitPoint();

    const scalar ppDist = pointPairDistance(edgePt);

    const vectorField& feNormals = feMesh.normals();
    const labelList& edNormalIs = feMesh.edgeNormals()[edHit.index()];

    if (edNormalIs.size() == 1)
    {
//        Info<< "Inserting open edge point group around " << edgePt << endl;
//        Info<< "    ppDist    = " << ppDist << nl
//            << "    edNormals = " << edNormalIs
//            << endl;

        const vector& n = feNormals[edNormalIs[0]];

        vector s = ppDist*(feMesh.edgeDirections()[edHit.index()] ^ n);

        plane facePlane(edgePt, n);

        Foam::point pt1 = edgePt + s + ppDist*n;
        Foam::point pt2 = edgePt - s + ppDist*n;

        Foam::point pt3 = reflectPointInPlane(pt1, facePlane);
        Foam::point pt4 = reflectPointInPlane(pt2, facePlane);

        pts.append(Vb(pt1, Vb::vtInternalSurface));
        pts.append(Vb(pt2, Vb::vtInternalSurface));
        pts.append(Vb(pt3, Vb::vtInternalSurface));
        pts.append(Vb(pt4, Vb::vtInternalSurface));
    }
    else
    {
        Info<< "NOT INSERTING OPEN EDGE POINT GROUP WITH MORE THAN 1 "
            << "EDGE NORMAL, NOT IMPLEMENTED" << endl;
    }
}


void Foam::conformalVoronoiMesh::createMultipleEdgePointGroup
(
    const extendedFeatureEdgeMesh& feMesh,
    const pointIndexHit& edHit,
    DynamicList<Vb>& pts
)
{
//    Info<< "NOT INSERTING MULTIPLE EDGE POINT GROUP, NOT IMPLEMENTED" << endl;

    const Foam::point& edgePt = edHit.hitPoint();

    const scalar ppDist = pointPairDistance(edgePt);

    const vectorField& feNormals = feMesh.normals();
    const labelList& edNormalIs = feMesh.edgeNormals()[edHit.index()];

    Info<< edNormalIs.size() << endl;

    // As this is a flat edge, there are two normals by definition
    const vector& nA = feNormals[edNormalIs[0]];
    const vector& nB = feNormals[edNormalIs[1]];

    // Average normal to remove any bias to one side, although as this
    // is a flat edge, the normals should be essentially the same
    const vector n = 0.5*(nA + nB);

    // Direction along the surface to the control point, sense of edge
    // direction not important, as +s and -s can be used because this
    // is a flat edge
    vector s = ppDist*(feMesh.edgeDirections()[edHit.index()] ^ n);

    createPointPair(ppDist, edgePt + s, n, pts);
    createPointPair(ppDist, edgePt - s, n, pts);
}


void Foam::conformalVoronoiMesh::reinsertFeaturePoints(bool distribute)
{
    Info<< nl << "Reinserting stored feature points" << endl;

    label preReinsertionSize(number_of_vertices());

    insertPoints(featureVertices_, distribute);

    const label nReinserted = returnReduce
    (
        label(number_of_vertices()) - preReinsertionSize,
        sumOp<label>()
    );

    Info<< "    Reinserted " << nReinserted << " vertices" << endl;
}


void Foam::conformalVoronoiMesh::createMixedFeaturePoints
(
    DynamicList<Vb>& pts
)
{
    if (foamyHexMeshControls().mixedFeaturePointPPDistanceCoeff() < 0)
    {
        Info<< nl << "Skipping specialised handling for mixed feature points"
            << endl;
        return;
    }

    const PtrList<extendedFeatureEdgeMesh>& feMeshes
    (
        geometryToConformTo_.features()
    );

    forAll(feMeshes, i)
    {
        const extendedFeatureEdgeMesh& feMesh = feMeshes[i];
        const labelListList& pointsEdges = feMesh.pointEdges();
        const pointField& points = feMesh.points();

        for
        (
            label ptI = feMesh.mixedStart();
            ptI < feMesh.nonFeatureStart();
            ptI++
        )
        {
            const Foam::point& featPt = points[ptI];

            if (!positionOnThisProc(featPt))
            {
                continue;
            }

            const labelList& pEds = pointsEdges[ptI];

            pointFeatureEdgesTypes pFEdgeTypes(ptI);

            const List<extendedFeatureEdgeMesh::edgeStatus> allEdStat
                = calcPointFeatureEdgesTypes(feMesh, pEds, pFEdgeTypes);

            bool specialisedSuccess = false;

            if (foamyHexMeshControls().specialiseFeaturePoints())
            {
                specialisedSuccess = createSpecialisedFeaturePoint
                (
                    feMesh, pEds, pFEdgeTypes, allEdStat, ptI, pts
                );
            }

            if (!specialisedSuccess && foamyHexMeshControls().edgeAiming())
            {
                // Specialisations available for some mixed feature points.  For
                // non-specialised feature points, inserting mixed internal and
                // external edge groups at feature point.

                // Skipping unsupported mixed feature point types
                bool skipEdge = false;

//                forAll(pEds, e)
//                {
//                    const label edgeI = pEds[e];
//
//                    const extendedFeatureEdgeMesh::edgeStatus edStatus
//                        = feMesh.getEdgeStatus(edgeI);
//
//                    if
//                    (
//                        edStatus == extendedFeatureEdgeMesh::OPEN
//                     || edStatus == extendedFeatureEdgeMesh::MULTIPLE
//                    )
//                    {
//                        Info<< "Edge type " << edStatus
//                            << " found for mixed feature point " << ptI
//                            << ". Not supported."
//                            << endl;
//
//                        skipEdge = true;
//                    }
//                }
//
//                if (skipEdge)
//                {
//                    Info<< "Skipping point " << ptI << nl << endl;
//
//                    continue;
//                }

//                createFeaturePoints(feMesh, ptI, pts, types);

                const Foam::point pt = points[ptI];

                const scalar edgeGroupDistance = mixedFeaturePointDistance(pt);

                forAll(pEds, e)
                {
                    const label edgeI = pEds[e];

                    const Foam::point edgePt =
                        pt + edgeGroupDistance*feMesh.edgeDirection(edgeI, ptI);

                    const pointIndexHit edgeHit(true, edgePt, edgeI);

                    createEdgePointGroup(feMesh, edgeHit, pts);
                }
            }
        }
    }
}
//
//
//void Foam::conformalVoronoiMesh::createFeaturePoints
//(
//    DynamicList<Foam::point>& pts,
//    DynamicList<label>& indices,
//    DynamicList<label>& types
//)
//{
//    const PtrList<extendedFeatureEdgeMesh>& feMeshes
//    (
//        geometryToConformTo_.features()
//    );
//
//    forAll(feMeshes, i)
//    {
//        const extendedFeatureEdgeMesh& feMesh(feMeshes[i]);
//
//        for
//        (
//            label ptI = feMesh.convexStart();
//            ptI < feMesh.mixedStart();
//            ++ptI
//        )
//        {
//            const Foam::point& featPt = feMesh.points()[ptI];
//
//            if (!positionOnThisProc(featPt))
//            {
//                continue;
//            }
//
//            const scalar searchRadiusSqr = 5.0*targetCellSize(featPt);
//
//            labelList indices = surfacePtLocationTreePtr_().findSphere
//            (
//                featPt,
//                searchRadiusSqr
//            );
//
//            pointField nearestSurfacePoints(indices.size());
//
//            forAll(indices, pI)
//            {
//                nearestSurfacePoints[pI] =
//                    surfaceConformationVertices_[indices[pI]];
//            }
//
//            forAll(feMesh.)
//
//            // Now find the nearest points within the edge cones.
//
//            // Calculate preliminary surface point locations
//
//
//        }
//    }
//}


void Foam::conformalVoronoiMesh::insertFeaturePoints()
{
    Info<< nl << "Conforming to feature points" << endl;

    Info<< "    Circulating edges is: "
        << foamyHexMeshControls().circulateEdges().asText() << nl
        << "    Guarding feature points is: "
        << foamyHexMeshControls().guardFeaturePoints().asText() << nl
        << "    Snapping to feature points is: "
        << foamyHexMeshControls().snapFeaturePoints().asText() << nl
        << "    Specialising feature points is: "
        << foamyHexMeshControls().specialiseFeaturePoints().asText() << endl;


    DynamicList<Vb> pts;

    const label preFeaturePointSize = number_of_vertices();

    createFeaturePoints(pts);

    createMixedFeaturePoints(pts);

    // Points added using the createEdgePointGroup function will be labelled as
    // internal/external feature edges. Relabel them as feature points,
    // otherwise they are inserted as both feature points and surface points.
    forAll(pts, pI)
    {
        Vb& pt = pts[pI];

        //if (pt.featureEdgePoint())
        {
            if (pt.internalBoundaryPoint())
            {
                pt.type() = Vb::vtInternalFeaturePoint;
            }
            else if (pt.externalBoundaryPoint())
            {
                pt.type() = Vb::vtExternalFeaturePoint;
            }
        }
    }

    // Insert the created points, distributing to the appropriate processor
    insertPoints(pts, true);

    if (foamyHexMeshControls().objOutput())
    {
        writePoints("featureVertices.obj", pts);
    }

    label nFeatureVertices = number_of_vertices() - preFeaturePointSize;
    reduce(nFeatureVertices, sumOp<label>());

    Info<< "    Inserted " << nFeatureVertices << " feature vertices" << endl;

    featureVertices_.clear();
    featureVertices_.setSize(pts.size());

    forAll(pts, pI)
    {
        featureVertices_[pI] = pts[pI];
    }

    constructFeaturePointLocations();
}


void Foam::conformalVoronoiMesh::constructFeaturePointLocations()
{
    DynamicList<Foam::point> ftPtLocs;

    const PtrList<extendedFeatureEdgeMesh>& feMeshes
    (
        geometryToConformTo_.features()
    );

    forAll(feMeshes, i)
    {
        const extendedFeatureEdgeMesh& feMesh(feMeshes[i]);

        if (foamyHexMeshControls().mixedFeaturePointPPDistanceCoeff() < 0)
        {
            // Ignoring mixed feature points
            for
            (
                label ptI = feMesh.convexStart();
                ptI < feMesh.mixedStart();
                ptI++
            )
            {
                ftPtLocs.append(feMesh.points()[ptI]);
            }
        }
        else
        {
            for
            (
                label ptI = feMesh.convexStart();
                ptI < feMesh.nonFeatureStart();
                ptI++
            )
            {
                ftPtLocs.append(feMesh.points()[ptI]);
            }
        }
    }

    featurePointLocations_.transfer(ftPtLocs);
}


Foam::List<Foam::pointIndexHit>
Foam::conformalVoronoiMesh::findSurfacePtLocationsNearFeaturePoint
(
    const Foam::point& featurePoint
) const
{
    DynamicList<pointIndexHit> dynPointList;

    const scalar searchRadiusSqr = 3*targetCellSize(featurePoint);

    labelList surfacePtList = surfacePtLocationTreePtr_().findSphere
    (
       featurePoint,
       searchRadiusSqr
    );

    forAll(surfacePtList, elemI)
    {
       label index = surfacePtList[elemI];

       const Foam::point& p
           = surfacePtLocationTreePtr_().shapes().shapePoints()[index];

       pointIndexHit nearHit(true, p, index);

       dynPointList.append(nearHit);
    }

    return dynPointList.shrink();
}


void Foam::conformalVoronoiMesh::addMasterAndSlavePoints
(
    const DynamicList<Foam::point>& masterPoints,
    const DynamicList<Foam::indexedVertexEnum::vertexType>& masterPointsTypes,
    const Map<DynamicList<autoPtr<plane> > >& masterPointReflections,
    DynamicList<Vb>& pts,
    const label ptI
) const
{
    typedef DynamicList<autoPtr<plane> >        planeDynList;
    typedef Foam::indexedVertexEnum::vertexType vertexType;

    forAll(masterPoints, pI)
    {
        // Append master to the list of points

        const Foam::point& masterPt = masterPoints[pI];
        const vertexType masterType = masterPointsTypes[pI];

//        Info<< "    Master = " << masterPt << endl;

        pts.append
        (
            Vb
            (
                masterPt,
                vertexCount() + pts.size(),
                masterType,
                Pstream::myProcNo()
            )
        );

        const label masterIndex = pts[pts.size() - 1].index();

        //meshTools::writeOBJ(strMasters, masterPt);

        const planeDynList& masterPointPlanes = masterPointReflections[pI];

        forAll(masterPointPlanes, planeI)
        {
            // Reflect master points in the planes and insert the slave points

            const plane& reflPlane = masterPointPlanes[planeI]();

            const Foam::point slavePt =
                reflectPointInPlane(masterPt, reflPlane);

//            Info<< "        Slave " << planeI << " = " << slavePt << endl;

            const vertexType slaveType =
            (
                masterType == Vb::vtInternalFeaturePoint
              ? Vb::vtExternalFeaturePoint // true
              : Vb::vtInternalFeaturePoint // false
            );

            pts.append
            (
                Vb
                (
                    slavePt,
                    vertexCount() + pts.size(),
                    slaveType,
                    Pstream::myProcNo()
                )
            );


            //meshTools::writeOBJ(strSlaves, slavePt);
        }
    }
}


Foam::label Foam::conformalVoronoiMesh::getSign
(
    const extendedFeatureEdgeMesh::edgeStatus eStatus
) const
{
    if (eStatus == extendedFeatureEdgeMesh::EXTERNAL)
    {
        return -1;
    }
    else if (eStatus == extendedFeatureEdgeMesh::INTERNAL)
    {
        return 1;
    }

    return 0;
}


void Foam::conformalVoronoiMesh::createMasterAndSlavePoints
(
    const extendedFeatureEdgeMesh& feMesh,
    const label ptI,
    DynamicList<Vb>& pts
) const
{
    typedef DynamicList<autoPtr<plane> >                planeDynList;
    typedef Foam::indexedVertexEnum::vertexType         vertexType;
    typedef Foam::extendedFeatureEdgeMesh::edgeStatus   edgeStatus;

    const Foam::point& featPt = feMesh.points()[ptI];

    if (!positionOnThisProc(featPt) || geometryToConformTo_.outside(featPt))
    {
        return;
    }

    const scalar ppDist = pointPairDistance(featPt);

    // Maintain a list of master points and the planes to relect them in
    DynamicList<Foam::point> masterPoints;
    DynamicList<vertexType> masterPointsTypes;
    Map<planeDynList> masterPointReflections;

    const labelList& featPtEdges = feMesh.featurePointEdges()[ptI];

    pointFeatureEdgesTypes pointEdgeTypes(ptI);

    const List<extendedFeatureEdgeMesh::edgeStatus> allEdStat
        = calcPointFeatureEdgesTypes(feMesh, featPtEdges, pointEdgeTypes);

//    Info<< nl << featPt << "  " << pointEdgeTypes;

    const_circulator<labelList> circ(featPtEdges);

    // Loop around the edges of the feature point
    if (circ.size()) do
    {
        const edgeStatus eStatusPrev = feMesh.getEdgeStatus(circ.prev());
        const edgeStatus eStatusCurr = feMesh.getEdgeStatus(circ());
//        const edgeStatus eStatusNext = feMesh.getEdgeStatus(circ.next());

//        Info<< "    Prev = "
//            << extendedFeatureEdgeMesh::edgeStatusNames_[eStatusPrev]
//            << " Curr = "
//            << extendedFeatureEdgeMesh::edgeStatusNames_[eStatusCurr]
////            << " Next = "
////            << extendedFeatureEdgeMesh::edgeStatusNames_[eStatusNext]
//            << endl;

        // Get the direction in which to move the point in relation to the
        // feature point
        label sign = getSign(eStatusCurr);

        const vector n = sharedFaceNormal(feMesh, circ(), circ.next());

        const vector pointMotionDirection = sign*0.5*ppDist*n;

//        Info<< "    Shared face normal      = " << n << endl;
//        Info<< "    Direction to move point = " << pointMotionDirection
//            << endl;

        if (masterPoints.empty())
        {
            // Initialise with the first master point
            Foam::point pt = featPt + pointMotionDirection;

            planeDynList firstPlane;
            firstPlane.append(autoPtr<plane>(new plane(featPt, n)));

            masterPoints.append(pt);

            masterPointsTypes.append
            (
                sign == 1
              ? Vb::vtExternalFeaturePoint // true
              : Vb::vtInternalFeaturePoint // false
            );

            //Info<< "    " << " " << firstPlane << endl;

//            const Foam::point reflectedPoint = reflectPointInPlane
//            (
//                masterPoints.last(),
//                firstPlane.last()()
//            );

            masterPointReflections.insert
            (
                masterPoints.size() - 1,
                firstPlane
            );
        }
//        else if
//        (
//            eStatusPrev == extendedFeatureEdgeMesh::INTERNAL
//         && eStatusCurr == extendedFeatureEdgeMesh::EXTERNAL
//        )
//        {
//            // Insert a new master point.
//            Foam::point pt = featPt + pointMotionDirection;
//
//            planeDynList firstPlane;
//            firstPlane.append(autoPtr<plane>(new plane(featPt, n)));
//
//            masterPoints.append(pt);
//
//            masterPointsTypes.append
//            (
//                sign == 1
//              ? Vb::vtExternalFeaturePoint // true
//              : Vb::vtInternalFeaturePoint // false
//            );
//
//            masterPointReflections.insert
//            (
//                masterPoints.size() - 1,
//                firstPlane
//            );
//        }
//        else if
//        (
//            eStatusPrev == extendedFeatureEdgeMesh::EXTERNAL
//         && eStatusCurr == extendedFeatureEdgeMesh::INTERNAL
//        )
//        {
//
//        }
        else
        {
            // Just add this face contribution to the latest master point

            masterPoints.last() += pointMotionDirection;

            masterPointReflections[masterPoints.size() - 1].append
            (
                autoPtr<plane>(new plane(featPt, n))
            );
        }

    } while (circ.circulate(CirculatorBase::CLOCKWISE));

    addMasterAndSlavePoints
    (
        masterPoints,
        masterPointsTypes,
        masterPointReflections,
        pts,
        ptI
    );
}


void Foam::conformalVoronoiMesh::createFeaturePoints(DynamicList<Vb>& pts)
{
    const PtrList<extendedFeatureEdgeMesh>& feMeshes
    (
        geometryToConformTo_.features()
    );

    forAll(feMeshes, i)
    {
        Info<< indent << "Edge mesh = " << feMeshes[i].name() << nl << endl;

        const extendedFeatureEdgeMesh& feMesh(feMeshes[i]);

        for
        (
            label ptI = feMesh.convexStart();
            ptI < feMesh.mixedStart();
//            ptI < feMesh.nonFeatureStart();
            ptI++
        )
        {
            createMasterAndSlavePoints(feMesh, ptI, pts);
        }

        if (foamyHexMeshControls().guardFeaturePoints())
        {
            for
            (
                //label ptI = feMesh.convexStart();
                label ptI = feMesh.mixedStart();
                ptI < feMesh.nonFeatureStart();
                ptI++
            )
            {
                pts.append
                (
                    Vb
                    (
                        feMesh.points()[ptI],
                        Vb::vtConstrained
                    )
                );
            }
        }
    }
}


//Foam::scalar Foam::conformalVoronoiMesh::pyramidVolume
//(
//    const Foam::point& apex,
//    const Foam::point& a,
//    const Foam::point& b,
//    const Foam::point& c,
//    const bool printInfo
//) const
//{
//    triPointRef tri(a, b, c);
//
//    tetPointRef tet(tri.a(), tri.b(), tri.c(), apex);
//
//    scalar volume = tet.mag();
//
////    scalar volume = (1.0/3.0)*constant::mathematical::pi;
////
////    K::Circle_3 circle(toPoint(a), toPoint(b), toPoint(c));
////
////    scalar height = mag(topoint(circle.center()) - apex);
////
////    volume *= circle.squared_radius()*height;
//
//    if (printInfo)
//    {
//        Info<< "Calculating volume of pyramid..." << nl
//            << "    Apex      : " << apex << nl
//            << "    Point a   : " << a << nl
//            << "    Point b   : " << b << nl
//            << "    Point c   : " << c << nl
//            << "        Center  : " << tri.centre() << nl
//            << "        Volume  : " << volume << endl;
//    }
//
//    return volume;
//}


//void Foam::conformalVoronoiMesh::createPyramidMasterPoint
//(
//    const Foam::point& apex,
//    const vectorField& edgeDirections,
//    Foam::point& masterPoint,
//    vectorField& norms
//) const
//{
//    pointField basePoints(edgeDirections.size() + 1);
//
//    forAll(edgeDirections, eI)
//    {
//        basePoints[eI] = edgeDirections[eI] + apex;
//    }
//
//    basePoints[edgeDirections.size() + 1] = apex;
//
//    face f(identity(edgeDirections.size()));
//
//    pyramidPointFaceRef p(f, apex);
//
//    const scalar ppDist = pointPairDistance(apex);
//
//
//    vector unitDir = f.centre();
//    unitDir /= mag(unitDir);
//
//    masterPoint = apex + ppDist*unitDir;
//
//    norms.setSize(edgeDirections.size());
//
//    forAll(norms, nI)
//    {
//        norms[nI] =
//    }
//}


//void Foam::conformalVoronoiMesh::createConvexConcaveFeaturePoints
//(
//    DynamicList<Foam::point>& pts,
//    DynamicList<label>& indices,
//    DynamicList<label>& types
//)
//{
//    const PtrList<extendedFeatureEdgeMesh>& feMeshes
//    (
//        geometryToConformTo_.features()
//    );
//
//    forAll(feMeshes, i)
//    {
//        const extendedFeatureEdgeMesh& feMesh(feMeshes[i]);
//
//        for
//        (
//            label ptI = feMesh.convexStart();
//            ptI < feMesh.mixedStart();
//            ptI++
//        )
//        {
//            const Foam::point& apex = feMesh.points()[ptI];
//
//            if (!positionOnThisProc(apex))
//            {
//                continue;
//            }
//
//            const vectorField& featPtEdgeDirections
//                = feMesh.featurePointEdgeDirections(ptI);
//
//            Foam::point masterPoint;
//            vectorField tetNorms;
//
//            createPyramidMasterPoint
//            (
//                apex,
//                featPtEdgeDirections,
//                masterPoint,
//                tetNorms
//            );
//
//
//
//            // Result when the points are eventually inserted (example n = 4)
//            // Add number_of_vertices() at insertion of first vertex to all
//            // numbers:
//            // pt           index type
//            // internalPt   0     1
//            // externalPt0  1     0
//            // externalPt1  2     0
//            // externalPt2  3     0
//            // externalPt3  4     0
//
//            // Result when the points are eventually inserted (example n = 5)
//            // Add number_of_vertices() at insertion of first vertex to all
//            // numbers:
//            // pt           index type
//            // internalPt0  0     5
//            // internalPt1  1     5
//            // internalPt2  2     5
//            // internalPt3  3     5
//            // internalPt4  4     5
//            // externalPt   5     4
//
//            if (geometryToConformTo_.inside(masterPoint))
//            {
//
//            }
//            else
//            {
//
//            }
//
//            pts.append(masterPoint);
//            indices.append(0);
//            types.append(1);
//
//            label internalPtIndex = -1;
//
//            forAll(tetNorms, nI)
//            {
//                const vector& n = tetNorms[nI];
//
//                Foam::point reflectedPoint
//                    = reflectPoint(featPt, masterPoint, n);
//
//                pts.append(reflectedPoint);
//                indices.append(0);
//                types.append(internalPtIndex--);
//            }
//        }
//    }
//}


Foam::vector Foam::conformalVoronoiMesh::sharedFaceNormal
(
    const extendedFeatureEdgeMesh& feMesh,
    const label edgeI,
    const label nextEdgeI
) const
{
    const labelList& edgeInormals = feMesh.edgeNormals()[edgeI];
    const labelList& nextEdgeInormals = feMesh.edgeNormals()[nextEdgeI];

    const vector& A1 = feMesh.normals()[edgeInormals[0]];
    const vector& A2 = feMesh.normals()[edgeInormals[1]];

    const vector& B1 = feMesh.normals()[nextEdgeInormals[0]];
    const vector& B2 = feMesh.normals()[nextEdgeInormals[1]];

//    Info<< "        A1 = " << A1 << endl;
//    Info<< "        A2 = " << A2 << endl;
//    Info<< "        B1 = " << B1 << endl;
//    Info<< "        B2 = " << B2 << endl;

    const scalar A1B1 = mag((A1 & B1) - 1.0);
    const scalar A1B2 = mag((A1 & B2) - 1.0);
    const scalar A2B1 = mag((A2 & B1) - 1.0);
    const scalar A2B2 = mag((A2 & B2) - 1.0);

//    Info<< "      A1B1 = " << A1B1 << endl;
//    Info<< "      A1B2 = " << A1B2 << endl;
//    Info<< "      A2B1 = " << A2B1 << endl;
//    Info<< "      A2B2 = " << A2B2 << endl;

    if (A1B1 < A1B2 && A1B1 < A2B1 && A1B1 < A2B2)
    {
        return 0.5*(A1 + B1);
    }
    else if (A1B2 < A1B1 && A1B2 < A2B1 && A1B2 < A2B2)
    {
        return 0.5*(A1 + B2);
    }
    else if (A2B1 < A1B1 && A2B1 < A1B2 && A2B1 < A2B2)
    {
        return 0.5*(A2 + B1);
    }
    else
    {
        return 0.5*(A2 + B2);
    }
}


Foam::List<Foam::point> Foam::conformalVoronoiMesh::reflectPointInPlanes
(
    const Foam::point p,
    const DynamicList<autoPtr<plane> >& planes
) const
{
    List<Foam::point> reflectedPoints(planes.size());

    forAll(planes, planeI)
    {
        reflectedPoints[planeI] = reflectPointInPlane(p, planes[planeI]());
    }

    return reflectedPoints;
}


Foam::point Foam::conformalVoronoiMesh::reflectPointInPlane
(
    const Foam::point p,
    const plane& planeN
) const
{
    const vector reflectedPtDir = p - planeN.nearestPoint(p);

    if ((planeN.normal() & reflectedPtDir) > 0)
    {
        return p - 2.0*planeN.distance(p)*planeN.normal();
    }
    else
    {
        return p + 2.0*planeN.distance(p)*planeN.normal();
    }
}
