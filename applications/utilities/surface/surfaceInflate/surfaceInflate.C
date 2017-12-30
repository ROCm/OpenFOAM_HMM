/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

Application
    surfaceInflate

Group
    grpSurfaceUtilities

Description
    Inflates surface. WIP. Checks for overlaps and locally lowers
    inflation distance

Usage
    - surfaceInflate [OPTION]

    \param -checkSelfIntersection \n
    Includes checks for self-intersection.

    \param -nSmooth
    Specify number of smoothing iterations

    \param -featureAngle
    Specify a feature angle


    E.g. inflate surface by 20mm with 1.5 safety factor:
        surfaceInflate DTC-scaled.obj 0.02 1.5 -featureAngle 45 -nSmooth 2

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "triSurfaceFields.H"
#include "triSurfaceMesh.H"
#include "triSurfaceGeoMesh.H"
#include "PackedBoolList.H"
#include "OBJstream.H"
#include "surfaceFeatures.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar calcVertexNormalWeight
(
    const triFace& f,
    const label pI,
    const vector& fN,
    const pointField& points
)
{
    label index = f.find(pI);

    if (index == -1)
    {
        FatalErrorInFunction
            << "Point not in face" << abort(FatalError);
    }

    const vector e1 = points[f[index]] - points[f[f.fcIndex(index)]];
    const vector e2 = points[f[index]] - points[f[f.rcIndex(index)]];

    return mag(fN)/(magSqr(e1)*magSqr(e2) + VSMALL);
}


tmp<vectorField> calcVertexNormals(const triSurface& surf)
{
    // Weighted average of normals of faces attached to the vertex
    // Weight = fA / (mag(e0)^2 * mag(e1)^2);
    tmp<vectorField> tpointNormals
    (
        new pointField(surf.nPoints(), Zero)
    );
    vectorField& pointNormals = tpointNormals.ref();

    const pointField& points = surf.points();
    const labelListList& pointFaces = surf.pointFaces();
    const labelList& meshPoints = surf.meshPoints();

    forAll(pointFaces, pI)
    {
        const labelList& pFaces = pointFaces[pI];

        forAll(pFaces, fI)
        {
            const label faceI = pFaces[fI];
            const triFace& f = surf[faceI];

            vector fN = f.normal(points);

            scalar weight = calcVertexNormalWeight
            (
                f,
                meshPoints[pI],
                fN,
                points
            );

            pointNormals[pI] += weight*fN;
        }

        pointNormals[pI] /= mag(pointNormals[pI]) + VSMALL;
    }

    return tpointNormals;
}


tmp<vectorField> calcPointNormals
(
    const triSurface& s,
    const PackedBoolList& isFeaturePoint,
    const List<surfaceFeatures::edgeStatus>& edgeStat
)
{
    //const pointField pointNormals(s.pointNormals());
    tmp<vectorField> tpointNormals(calcVertexNormals(s));
    vectorField& pointNormals = tpointNormals.ref();


    // feature edges: create edge normals from edgeFaces only.
    {
        const labelListList& edgeFaces = s.edgeFaces();

        labelList nNormals(s.nPoints(), 0);
        forAll(edgeStat, edgeI)
        {
            if (edgeStat[edgeI] != surfaceFeatures::NONE)
            {
                const edge& e = s.edges()[edgeI];
                forAll(e, i)
                {
                    if (!isFeaturePoint[e[i]])
                    {
                        pointNormals[e[i]] = Zero;
                    }
                }
            }
        }

        forAll(edgeStat, edgeI)
        {
            if (edgeStat[edgeI] != surfaceFeatures::NONE)
            {
                const labelList& eFaces = edgeFaces[edgeI];

                // Get average edge normal
                vector n = Zero;
                forAll(eFaces, i)
                {
                    n += s.faceNormals()[eFaces[i]];
                }
                n /= eFaces.size();


                // Sum to points
                const edge& e = s.edges()[edgeI];
                forAll(e, i)
                {
                    if (!isFeaturePoint[e[i]])
                    {
                        pointNormals[e[i]] += n;
                        nNormals[e[i]]++;
                    }
                }
            }
        }

        forAll(nNormals, pointI)
        {
            if (nNormals[pointI] > 0)
            {
                pointNormals[pointI] /= mag(pointNormals[pointI]);
            }
        }
    }


    forAll(pointNormals, pointI)
    {
        if (mag(mag(pointNormals[pointI])-1) > SMALL)
        {
            FatalErrorInFunction
                << "unitialised"
                << exit(FatalError);
        }
    }

    return tpointNormals;
}


void detectSelfIntersections
(
    const triSurfaceMesh& s,
    PackedBoolList& isEdgeIntersecting
)
{
    const edgeList& edges = s.edges();
    const indexedOctree<treeDataTriSurface>& tree = s.tree();
    const labelList& meshPoints = s.meshPoints();
    const pointField& points = s.points();

    isEdgeIntersecting.setSize(edges.size());
    isEdgeIntersecting = false;

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        pointIndexHit hitInfo
        (
            tree.findLine
            (
                points[meshPoints[e[0]]],
                points[meshPoints[e[1]]],
                treeDataTriSurface::findSelfIntersectOp(tree, edgeI)
            )
        );

        if (hitInfo.hit())
        {
            isEdgeIntersecting[edgeI] = true;
        }
    }
}


label detectIntersectionPoints
(
    const scalar tolerance,
    const triSurfaceMesh& s,

    const scalar extendFactor,
    const pointField& initialPoints,
    const vectorField& displacement,

    const bool checkSelfIntersect,
    const PackedBoolList& initialIsEdgeIntersecting,

    PackedBoolList& isPointOnHitEdge,
    scalarField& scale
)
{
    const pointField initialLocalPoints(initialPoints, s.meshPoints());
    const List<labelledTri>& localFaces = s.localFaces();


    label nHits = 0;
    isPointOnHitEdge.setSize(s.nPoints());
    isPointOnHitEdge = false;


    // 1. Extrusion offset vectors intersecting new surface location
    {
        scalar tol = max(tolerance, 10*s.tolerance());

        // Collect all the edge vectors. Slightly shorten the edges to prevent
        // finding lots of intersections. The fast triangle intersection routine
        // has problems with rays passing through a point of the
        // triangle.
        pointField start(initialLocalPoints+tol*displacement);
        pointField end(initialLocalPoints+extendFactor*displacement);

        List<pointIndexHit> hits;
        s.findLineAny(start, end, hits);

        forAll(hits, pointI)
        {
            if
            (
                hits[pointI].hit()
            &&  !localFaces[hits[pointI].index()].found(pointI)
            )
            {
                scale[pointI] = max(0.0, scale[pointI]-0.2);

                isPointOnHitEdge[pointI] = true;
                nHits++;
            }
        }
    }


    // 2. (new) surface self intersections
    if (checkSelfIntersect)
    {
        PackedBoolList isEdgeIntersecting;
        detectSelfIntersections(s, isEdgeIntersecting);

        const edgeList& edges = s.edges();
        const pointField& points = s.points();

        forAll(edges, edgeI)
        {
            const edge& e = edges[edgeI];

            if (isEdgeIntersecting[edgeI] && !initialIsEdgeIntersecting[edgeI])
            {
                if (isPointOnHitEdge.set(e[0]))
                {
                    label start = s.meshPoints()[e[0]];
                    const point& pt = points[start];

                    Pout<< "Additional self intersection at "
                        << pt
                        << endl;

                    scale[e[0]] = max(0.0, scale[e[0]]-0.2);
                    nHits++;
                }
                if (isPointOnHitEdge.set(e[1]))
                {
                    label end = s.meshPoints()[e[1]];
                    const point& pt = points[end];

                    Pout<< "Additional self intersection at "
                        << pt
                        << endl;

                    scale[e[1]] = max(0.0, scale[e[1]]-0.2);
                    nHits++;
                }
            }
        }
    }

    return nHits;
}


tmp<scalarField> avg
(
    const triSurface& s,
    const scalarField& fld,
    const scalarField& edgeWeights
)
{
    tmp<scalarField> tres(new scalarField(s.nPoints(), 0.0));
    scalarField& res = tres.ref();

    scalarField sumWeight(s.nPoints(), 0.0);

    const edgeList& edges = s.edges();

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];
        const scalar w = edgeWeights[edgeI];

        res[e[0]] += w*fld[e[1]];
        sumWeight[e[0]] += w;

        res[e[1]] += w*fld[e[0]];
        sumWeight[e[1]] += w;
    }

    // Average
    // ~~~~~~~

    forAll(res, pointI)
    {
        if (mag(sumWeight[pointI]) < VSMALL)
        {
            // Unconnected point. Take over original value
            res[pointI] = fld[pointI];
        }
        else
        {
            res[pointI] /= sumWeight[pointI];
        }
    }

    return tres;
}


void minSmooth
(
    const triSurface& s,
    const PackedBoolList& isAffectedPoint,
    const scalarField& fld,
    scalarField& newFld
)
{

    const edgeList& edges = s.edges();
    const pointField& points = s.points();
    const labelList& mp = s.meshPoints();

    scalarField edgeWeights(edges.size());
    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];
        scalar w = mag(points[mp[e[0]]]-points[mp[e[1]]]);

        edgeWeights[edgeI] = 1.0/(max(w, SMALL));
    }

    tmp<scalarField> tavgFld = avg(s, fld, edgeWeights);

    const scalarField& avgFld = tavgFld();

    forAll(fld, pointI)
    {
        if (isAffectedPoint.get(pointI) == 1)
        {
            newFld[pointI] = min
            (
                fld[pointI],
                0.5*fld[pointI] + 0.5*avgFld[pointI]
                //avgFld[pointI]
            );
        }
    }
}


void lloydsSmoothing
(
    const label nSmooth,
    triSurface& s,
    const PackedBoolList& isFeaturePoint,
    const List<surfaceFeatures::edgeStatus>& edgeStat,
    const PackedBoolList& isAffectedPoint
)
{
    const labelList& meshPoints = s.meshPoints();
    const edgeList& edges = s.edges();


    PackedBoolList isSmoothPoint(isAffectedPoint);
    // Extend isSmoothPoint
    {
        PackedBoolList newIsSmoothPoint(isSmoothPoint);
        forAll(edges, edgeI)
        {
            const edge& e = edges[edgeI];
            if (isSmoothPoint[e[0]])
            {
                newIsSmoothPoint[e[1]] = true;
            }
            if (isSmoothPoint[e[1]])
            {
                newIsSmoothPoint[e[0]] = true;
            }
        }
        Info<< "From points-to-smooth " << isSmoothPoint.count()
            << " to " << newIsSmoothPoint.count()
            << endl;
        isSmoothPoint.transfer(newIsSmoothPoint);
    }

    // Do some smoothing (Lloyds algorithm) around problematic points
    for (label i = 0; i < nSmooth; i++)
    {
        const labelListList& pointFaces = s.pointFaces();
        const pointField& faceCentres = s.faceCentres();

        pointField newPoints(s.points());
        forAll(isSmoothPoint, pointI)
        {
            if (isSmoothPoint[pointI] && !isFeaturePoint[pointI])
            {
                const labelList& pFaces = pointFaces[pointI];

                point avg(Zero);
                forAll(pFaces, pFaceI)
                {
                    avg += faceCentres[pFaces[pFaceI]];
                }
                avg /= pFaces.size();
                newPoints[meshPoints[pointI]] = avg;
            }
        }


        // Move points on feature edges only according to feature edges.

        const pointField& points = s.points();

        vectorField pointSum(s.nPoints(), Zero);
        labelList nPointSum(s.nPoints(), 0);

        forAll(edges, edgeI)
        {
            if (edgeStat[edgeI] != surfaceFeatures::NONE)
            {
                const edge& e = edges[edgeI];
                const point& start = points[meshPoints[e[0]]];
                const point& end = points[meshPoints[e[1]]];
                point eMid = 0.5*(start+end);
                pointSum[e[0]] += eMid;
                nPointSum[e[0]]++;
                pointSum[e[1]] += eMid;
                nPointSum[e[1]]++;
            }
        }

        forAll(pointSum, pointI)
        {
            if
            (
                isSmoothPoint[pointI]
             && isFeaturePoint[pointI]
             && nPointSum[pointI] >= 2
            )
            {
                newPoints[meshPoints[pointI]] =
                    pointSum[pointI]/nPointSum[pointI];
            }
        }


        s.movePoints(newPoints);


        // Extend isSmoothPoint
        {
            PackedBoolList newIsSmoothPoint(isSmoothPoint);
            forAll(edges, edgeI)
            {
                const edge& e = edges[edgeI];
                if (isSmoothPoint[e[0]])
                {
                    newIsSmoothPoint[e[1]] = true;
                }
                if (isSmoothPoint[e[1]])
                {
                    newIsSmoothPoint[e[0]] = true;
                }
            }
            Info<< "From points-to-smooth " << isSmoothPoint.count()
                << " to " << newIsSmoothPoint.count()
                << endl;
            isSmoothPoint.transfer(newIsSmoothPoint);
        }
    }
}



// Main program:

int main(int argc, char *argv[])
{
    argList::addNote("Inflates surface according to point normals.");

    argList::noParallel();
    argList::addNote
    (
        "Creates inflated version of surface using point normals."
        " Takes surface, distance to inflate and additional safety factor"
    );
    argList::addBoolOption
    (
        "checkSelfIntersection",
        "also check for self-intersection"
    );
    argList::addOption
    (
        "nSmooth",
        "integer",
        "number of smoothing iterations (default 20)"
    );
    argList::addOption
    (
        "featureAngle",
        "scalar",
        "feature angle"
    );
    argList::addBoolOption
    (
        "debug",
        "switch on additional debug information"
    );

    argList::addArgument("inputFile");
    argList::addArgument("distance");
    argList::addArgument("safety factor [1..]");

    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();

    const word inputName(args[1]);
    const scalar distance(args.argRead<scalar>(2));
    const scalar extendFactor(args.argRead<scalar>(3));
    const bool checkSelfIntersect = args.optionFound("checkSelfIntersection");
    const label nSmooth = args.optionLookupOrDefault("nSmooth", 10);
    const scalar featureAngle = args.optionLookupOrDefault<scalar>
    (
        "featureAngle",
        180
    );
    const bool debug = args.optionFound("debug");


    Info<< "Inflating surface " << inputName
        << " according to point normals" << nl
        << "    distance               : " << distance << nl
        << "    safety factor          : " << extendFactor << nl
        << "    self intersection test : " << checkSelfIntersect << nl
        << "    smoothing iterations   : " << nSmooth << nl
        << "    feature angle          : " << featureAngle << nl
        << endl;


    if (extendFactor < 1 || extendFactor > 10)
    {
        FatalErrorInFunction
            << "Illegal safety factor " << extendFactor
            << ". It is usually 1..2"
            << exit(FatalError);
    }



    // Load triSurfaceMesh
    triSurfaceMesh s
    (
        IOobject
        (
            inputName,                          // name
            runTime.constant(),                 // instance
            "triSurface",                       // local
            runTime,                            // registry
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    );

    // Mark features
    const surfaceFeatures features(s, 180.0-featureAngle);

    Info<< "Detected features:" << nl
        << "    feature points : " << features.featurePoints().size()
        << " out of " << s.nPoints() << nl
        << "    feature edges : " << features.featureEdges().size()
        << " out of " << s.nEdges() << nl
        << endl;

    PackedBoolList isFeaturePoint(s.nPoints());
    forAll(features.featurePoints(), i)
    {
        label pointI = features.featurePoints()[i];
        isFeaturePoint[pointI] = true;
    }
    const List<surfaceFeatures::edgeStatus> edgeStat(features.toStatus());




    const pointField initialPoints(s.points());


    // Construct scale
    Info<< "Constructing field scale\n" << endl;
    triSurfacePointScalarField scale
    (
        IOobject
        (
            "scale",                            // name
            runTime.timeName(),                 // instance
            s,                                  // registry
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        s,
        dimensionedScalar("scale", dimLength, 1.0)
    );


    // Construct unit normals

    Info<< "Calculating vertex normals\n" << endl;
    const pointField pointNormals
    (
        calcPointNormals
        (
            s,
            isFeaturePoint,
            edgeStat
        )
    );


    // Construct pointDisplacement
    Info<< "Constructing field pointDisplacement\n" << endl;
    triSurfacePointVectorField pointDisplacement
    (
        IOobject
        (
            "pointDisplacement",                // name
            runTime.timeName(),                 // instance
            s,                                  // registry
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        s,
        dimLength,
        distance*scale*pointNormals
    );


    const labelList& meshPoints = s.meshPoints();


    // Any point on any intersected edge in any of the iterations
    PackedBoolList isScaledPoint(s.nPoints());


    // Detect any self intersections on initial mesh
    PackedBoolList initialIsEdgeIntersecting;
    if (checkSelfIntersect)
    {
        detectSelfIntersections(s, initialIsEdgeIntersecting);
    }
    else
    {
        // Mark all edges as already self intersecting so avoid detecting any
        // new ones
        initialIsEdgeIntersecting.setSize(s.nEdges(), true);
    }


    // Inflate
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Move to new location
        pointField newPoints(initialPoints);
        forAll(meshPoints, i)
        {
            newPoints[meshPoints[i]] += pointDisplacement[i];
        }
        s.movePoints(newPoints);
        Info<< "Bounding box : " << s.bounds() << endl;


        // Work out scale from pointDisplacement
        forAll(scale, pointI)
        {
            if (s.pointFaces()[pointI].size())
            {
                scale[pointI] = mag(pointDisplacement[pointI])/distance;
            }
            else
            {
                scale[pointI] = 1.0;
            }
        }


        Info<< "Scale        : " << gAverage(scale) << endl;
        Info<< "Displacement : " << gAverage(pointDisplacement) << endl;



        // Detect any intersections and scale back
        PackedBoolList isAffectedPoint;
        label nIntersections = detectIntersectionPoints
        (
            1e-9,       // intersection tolerance
            s,
            extendFactor,
            initialPoints,
            pointDisplacement,

            checkSelfIntersect,
            initialIsEdgeIntersecting,

            isAffectedPoint,
            scale
        );
        Info<< "Detected " << nIntersections << " intersections" << nl
            << endl;

        if (nIntersections == 0)
        {
            runTime.write();
            break;
        }


        // Accumulate all affected points
        forAll(isAffectedPoint, pointI)
        {
            if (isAffectedPoint[pointI])
            {
                isScaledPoint[pointI] = true;
            }
        }

        // Smear out lowering of scale so any edges not found are
        // still included
        for (label i = 0; i < nSmooth; i++)
        {
            triSurfacePointScalarField oldScale(scale);
            oldScale.rename("oldScale");
            minSmooth
            (
                s,
                PackedBoolList(s.nPoints(), true),
                oldScale,
                scale
            );
        }


        // From scale update the pointDisplacement
        pointDisplacement *= distance*scale/mag(pointDisplacement);


        // Do some smoothing (Lloyds algorithm)
        lloydsSmoothing(nSmooth, s, isFeaturePoint, edgeStat, isAffectedPoint);


        // Update pointDisplacement
        const pointField& pts = s.points();
        forAll(meshPoints, i)
        {
            label meshPointI = meshPoints[i];
            pointDisplacement[i] = pts[meshPointI]-initialPoints[meshPointI];
        }


        // Write
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }


    Info<< "Detected overall intersecting points " << isScaledPoint.count()
        << " out of " << s.nPoints() << nl << endl;


    if (debug)
    {
        OBJstream str(runTime.path()/"isScaledPoint.obj");
        forAll(isScaledPoint, pointI)
        {
            if (isScaledPoint[pointI])
            {
                str.write(initialPoints[meshPoints[pointI]]);
            }
        }
    }


    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
