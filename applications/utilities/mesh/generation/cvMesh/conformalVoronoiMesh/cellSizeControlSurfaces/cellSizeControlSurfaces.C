/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "cellSizeControlSurfaces.H"
#include "conformalVoronoiMesh.H"
#include "cellSizeFunction.H"
#include "triSurfaceMesh.H"
#include "tetrahedron.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(cellSizeControlSurfaces, 0);

}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

bool Foam::cellSizeControlSurfaces::evalCellSizeFunctions
(
    const point& pt,
    scalar& minSize
) const
{
    bool anyFunctionFound = false;

    // // Regions requesting with the same priority take the average

    // scalar sizeAccumulator = 0;
    // scalar numberOfFunctions = 0;

    // label previousPriority = defaultPriority_;

    // if (cellSizeFunctions_.size())
    // {
    //     previousPriority =
    //         cellSizeFunctions_[cellSizeFunctions_.size() - 1].priority();

    //     forAll(cellSizeFunctions_, i)
    //     {
    //         const cellSizeFunction& cSF = cellSizeFunctions_[i];

    //         if (cSF.priority() < previousPriority && numberOfFunctions > 0)
    //         {
    //             return sizeAccumulator/numberOfFunctions;
    //         }

    //         scalar sizeI;

    //         if (cSF.cellSize(pt, sizeI))
    //         {
    //             anyFunctionFound = true;

    //             previousPriority = cSF.priority();

    //             sizeAccumulator += sizeI;
    //             numberOfFunctions++;
    //         }
    //     }
    // }

    // if (previousPriority == defaultPriority_ || numberOfFunctions == 0)
    // {
    //     sizeAccumulator += defaultCellSize_;
    //     numberOfFunctions++;
    // }

    // minSize = sizeAccumulator/numberOfFunctions;

    // return anyFunctionFound;

    // Regions requesting with the same priority take the smallest

    if (cellSizeFunctions_.size())
    {
        // Maintain priority of current hit. Initialise so it always goes
        // through at least once.
        label previousPriority = -1;

        forAll(cellSizeFunctions_, i)
        {
            const cellSizeFunction& cSF = cellSizeFunctions_[i];

            if (debug)
            {
                Info<< "size function "
                    << allGeometry_.names()[surfaces_[i]]
                    << " priority " << cSF.priority()
                    << endl;
            }

            if (cSF.priority() < previousPriority)
            {
                return minSize;
            }

            scalar sizeI;

            if (cSF.cellSize(pt, sizeI))
            {
                anyFunctionFound = true;

                if (cSF.priority() == previousPriority)
                {
                    if (sizeI < minSize)
                    {
                        minSize = sizeI;
                    }
                }
                else
                {
                    minSize = sizeI;
                }

                if (debug)
                {
                    Info<< "sizeI " << sizeI << " minSize " << minSize << endl;
                }

                previousPriority = cSF.priority();
            }
        }
    }

    return anyFunctionFound;
}


bool Foam::cellSizeControlSurfaces::checkCoplanarTet
(
    Cell_handle c,
    const scalar tol
) const
{
    plane triPlane
    (
        topoint(c->vertex(0)->point()),
        topoint(c->vertex(1)->point()),
        topoint(c->vertex(2)->point())
    );

    // Check if the four points are roughly coplanar. If they are then we
    // cannot calculate the circumcentre. Better test might be the volume
    // of the tet.
    if (triPlane.distance(topoint(c->vertex(3)->point())) < tol)
    {
        return true;
    }

    return false;
}


bool Foam::cellSizeControlSurfaces::checkClosePoints
(
    Cell_handle c,
    const scalar tol
) const
{
    for (label v = 0; v < 4; ++v)
    {
        for (label vA = v + 1; vA < 4; ++vA)
        {
            if
            (
                mag
                (
                    topoint(c->vertex(v)->point())
                  - topoint(c->vertex(vA)->point())
                )
              < tol
            )
            {
                return true;
            }
        }
    }

    return false;
}


Foam::label Foam::cellSizeControlSurfaces::refineTriangulation
(
    const scalar factor
)
{
    // Check the tets and insert new points if necessary
    List<Foam::point> corners(4);
    List<scalar> values(4);

    DynamicList<Foam::point> pointsToInsert(T_.number_of_vertices());
    DynamicList<scalar> valuesToInsert(T_.number_of_vertices());

    for
    (
        Delaunay::Finite_cells_iterator c = T_.finite_cells_begin();
        c != T_.finite_cells_end();
        ++c
    )
    {
        //const point& newPoint = tet.centre();
        //Cell_handle ch = c;

//        Info<< count++ << endl;
//        Info<< "    " << topoint(c->vertex(0)->point()) << nl
//            << "    " << topoint(c->vertex(1)->point()) << nl
//            << "    " << topoint(c->vertex(2)->point()) << nl
//            << "    " << topoint(c->vertex(3)->point()) << endl;

        scalar minDist = 1e-6*defaultCellSize_;

        if (checkClosePoints(c, minDist) || checkCoplanarTet(c, minDist))
        {
            continue;
        }

        const Point circumcenter = CGAL::circumcenter
        (
            c->vertex(0)->point(),
            c->vertex(1)->point(),
            c->vertex(2)->point(),
            c->vertex(3)->point()
        );

        pointFromPoint newPoint = topoint(circumcenter);

        if (geometryToConformTo_.outside(newPoint))
        {
            continue;
        }

        Cell_handle ch = T_.locate
        (
            Point(newPoint.x(), newPoint.y(), newPoint.z())
        );

        forAll(corners, pI)
        {
            corners[pI] = topoint(ch->vertex(pI)->point());
            values[pI] = ch->vertex(pI)->value();
        }

        tetPointRef tet(corners[0], corners[1], corners[2], corners[3]);

        scalarList bary;
        tet.barycentric(newPoint, bary);

        scalar interpolatedSize = 0;
        forAll(bary, pI)
        {
            interpolatedSize += bary[pI]*ch->vertex(pI)->value();
        }

        // Find largest gradient
        label maxGradCorner = -1;
        scalar maxGradient = 0.0;
        forAll(corners, pI)
        {
            const scalar distance = mag(newPoint - corners[pI]);
            const scalar diffSize = interpolatedSize - values[pI];

            const scalar gradient = diffSize/distance;

            if (gradient > maxGradient)
            {
                maxGradient = gradient;
                maxGradCorner = pI;
            }
        }

//        if (wallSize < 0.5*defaultCellSize_)
//        {
//        Info<< "Centre             : " << centre
//            << " (Default Size: " << defaultCellSize_ << ")" << nl
//            << "Interpolated Size  : " << interpolatedSize << nl
//            << "Distance from wall : " << distanceSize << nl
//            << "Wall size          : " << wallSize << nl
//            << "distanceGradient   : " << distanceGradient << nl
//            << "interpGradient     : " << interpolatedGradient << nl << endl;
//        }

        scalar minCellSize = 1e-6;
        scalar initialMaxGradient = 0;//0.2*factor;
        scalar initialMinGradient = 0;//0.01*(1.0/factor);
        scalar idealGradient = 0.2;



        // Reduce strong gradients
        if (maxGradient > initialMaxGradient)
        {
            const scalar distance2 = mag(newPoint - corners[maxGradCorner]);

            scalar newSize
                = values[maxGradCorner] + idealGradient*distance2;

            pointsToInsert.append(newPoint);
            valuesToInsert.append
            (
                max(min(newSize, defaultCellSize_), minCellSize)
            );
        }
        else if (maxGradient < -initialMaxGradient)
        {
            const scalar distance2 = mag(newPoint - corners[maxGradCorner]);

            scalar newSize
                = values[maxGradCorner] - idealGradient*distance2;

            pointsToInsert.append(newPoint);
            valuesToInsert.append
            (
                max(min(newSize, defaultCellSize_), minCellSize)
            );
        }

        // Increase small gradients
        if
        (
            maxGradient < initialMinGradient
         && maxGradient > 0
         && interpolatedSize < 0.5*defaultCellSize_
        )
        {
            const scalar distance2 = mag(newPoint - corners[maxGradCorner]);

            scalar newSize
                = values[maxGradCorner] + idealGradient*distance2;

            pointsToInsert.append(newPoint);
            valuesToInsert.append
            (
                max(min(newSize, defaultCellSize_), minCellSize)
            );
        }
        else if
        (
            maxGradient > -initialMinGradient
         && maxGradient < 0
         && interpolatedSize > 0.5*defaultCellSize_
        )
        {
            const scalar distance2 = mag(newPoint - corners[maxGradCorner]);

            scalar newSize
                = values[maxGradCorner] - idealGradient*distance2;

            pointsToInsert.append(newPoint);
            valuesToInsert.append
            (
                max(min(newSize, defaultCellSize_), minCellSize)
            );
        }
    }

    if (!pointsToInsert.empty())
    {
        Info<< "    Minimum Cell Size : " << min(valuesToInsert) << nl
            << "    Average Cell Size : " << average(valuesToInsert) << nl
            << "    Maximum Cell Size : " << max(valuesToInsert) << endl;

        forAll(pointsToInsert, pI)
        {
            Foam::point p = pointsToInsert[pI];

            Vertex_handle v = T_.insert(Point(p.x(), p.y(), p.z()));

            v->value(valuesToInsert[pI]);
        }
    }

    return pointsToInsert.size();
}


void Foam::cellSizeControlSurfaces::writeRefinementTriangulation()
{
    OFstream str("refinementTriangulation.obj");

    label count = 0;

    Info<< "Write refinementTriangulation" << endl;

    for
    (
        Delaunay::Finite_edges_iterator e = T_.finite_edges_begin();
        e != T_.finite_edges_end();
        ++e
    )
    {
        Cell_handle c = e->first;
        Vertex_handle vA = c->vertex(e->second);
        Vertex_handle vB = c->vertex(e->third);

        pointFromPoint p1 = topoint(vA->point());
        pointFromPoint p2 = topoint(vB->point());

        meshTools::writeOBJ(str, p1, p2, count);
    }


    OFstream strPoints("refinementObject.obj");

    Info<< "Write refinementObject" << endl;

    for
    (
        Delaunay::Finite_vertices_iterator v = T_.finite_vertices_begin();
        v != T_.finite_vertices_end();
        ++v
    )
    {
        pointFromPoint p = topoint(v->point());

        meshTools::writeOBJ
        (
            strPoints,
            p,
            point(p.x() + v->value(), p.y(), p.z())
        );
    }

//    OFstream strDual("refinementDualPoints.obj");
//
//    Info<< "Write refinementDualPoints" << endl;
//
//    for
//    (
//        Delaunay::Finite_cells_iterator c = T_.finite_cells_begin();
//        c != T_.finite_cells_end();
//        ++c
//    )
//    {
//        Point circumcenter = CGAL::circumcenter
//        (
//            c->vertex(0)->point(),
//            c->vertex(1)->point(),
//            c->vertex(2)->point(),
//            c->vertex(3)->point()
//        );
//
//        pointFromPoint p = topoint(circumcenter);
//
//        if (geometryToConformTo_.inside(p))
//        {
//            meshTools::writeOBJ
//            (
//                strDual,
//                p
//            );
//        }
//    }

    if (T_.is_valid())
    {
        Info<< "    Triangulation is valid" << endl;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::cellSizeControlSurfaces::writeRefinementTriangulation()"
        )   << "Triangulation is not valid"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellSizeControlSurfaces::cellSizeControlSurfaces
(
    const searchableSurfaces& allGeometry,
    const conformationSurfaces& geometryToConformTo,
    const dictionary& motionControlDict
)
:
    allGeometry_(allGeometry),
    geometryToConformTo_(geometryToConformTo),
    surfaces_(),
    cellSizeFunctions_(),
    defaultCellSize_(readScalar(motionControlDict.lookup("defaultCellSize"))),
    defaultPriority_
    (
        motionControlDict.lookupOrDefault<label>("defaultPriority", 0)
    )
{
    const dictionary& surfacesDict
    (
        motionControlDict.subDict("cellSizeControlGeometry")
    );

    Info<< nl << "Reading cellSizeControlGeometry" << endl;

    surfaces_.setSize(surfacesDict.size());

    cellSizeFunctions_.setSize(surfacesDict.size());

    labelList priorities(surfacesDict.size());

    label surfI = 0;

    forAllConstIter(dictionary, surfacesDict, iter)
    {
        const dictionary& surfaceSubDict
        (
            surfacesDict.subDict(iter().keyword())
        );

        // If the "surface" keyword is not found in the dictionary, assume that
        // the name of the dictionary is the surface.  Distinction required to
        // allow the same surface to be used multiple times to supply multiple
        // cellSizeFunctions

        word surfaceName = surfaceSubDict.lookupOrDefault<word>
        (
            "surface",
            iter().keyword()
        );

        surfaces_[surfI] = allGeometry_.findSurfaceID(surfaceName);

        if (surfaces_[surfI] < 0)
        {
            FatalErrorIn
            (
                "Foam::cellSizeControlSurfaces::cellSizeControlSurfaces"
            )   << "No surface " << surfaceName << " found. "
                << "Valid geometry is " << nl << allGeometry_.names()
                << exit(FatalError);
        }

        const searchableSurface& surface = allGeometry_[surfaces_[surfI]];

        Info<< nl << "    " << iter().keyword() << nl
            << "    surface: " << surfaceName << nl
            << "    size   : " << surface.size() << endl;

        cellSizeFunctions_.set
        (
            surfI,
            cellSizeFunction::New
            (
                surfaceSubDict,
                surface
            )
        );

        priorities[surfI] = cellSizeFunctions_[surfI].priority();

        surfI++;

        if (isA<triSurfaceMesh>(surface))
        {
            const triSurfaceMesh& tsm
                = refCast<const triSurfaceMesh>(surface);

            const pointField& points = tsm.points();

            Info<< "    number of points: " << tsm.nPoints() << endl;

            std::vector<std::pair<Point, scalar> > pointsToInsert;

            forAll(points, pI)
            {
                size_t nVert = T_.number_of_vertices();

                Vertex_handle v = T_.insert
                (
                    Point(points[pI].x(), points[pI].y(), points[pI].z())
                );

                if (T_.number_of_vertices() == nVert)
                {
                    Info<< "Failed to insert point : " << points[pI] << endl;
                }

                // Get the value of the point from surfaceCellSizeFunction. If
                // adding points internally then will need to interpolate.
                scalar newSize = 0;
                cellSizeFunctions_[surfI-1].cellSize(points[pI], newSize);
                v->value(newSize);
            }
        }
    }

    scalar factor = 1.0;
    label maxIteration = 1;

    for (label iteration = 0; iteration < maxIteration; ++iteration)
    {
        Info<< "Iteration : " << iteration << endl;

        label nRefined = refineTriangulation(factor);

        Info<< "    Number of cells refined in refinement iteration : "
            << nRefined << nl << endl;

        if (nRefined <= 0 && iteration != 0)
        {
            break;
        }

        factor *= 1.5;
    }

    writeRefinementTriangulation();

    Info<< nl << "Refinement triangulation information: " << endl;
    Info<< "    Number of vertices: " << label(T_.number_of_vertices()) << endl;
    Info<< "    Number of cells   : "
        << label(T_.number_of_finite_cells()) << endl;
    Info<< "    Number of faces   : "
        << label(T_.number_of_finite_facets()) << endl;
    Info<< "    Number of edges   : "
        << label(T_.number_of_finite_edges()) << endl;
    Info<< "    Dimensionality    : " << label(T_.dimension()) << nl << endl;


    // Sort cellSizeFunctions_ and surfaces_ by priority.  Cut off any surfaces
    // where priority < defaultPriority_


    labelList sortedIndices;

    sortedOrder(priorities, sortedIndices);

    sortedIndices = invert(sortedIndices.size(), sortedIndices);

    // Reverse the sort order
    sortedIndices = (sortedIndices.size() - 1) - sortedIndices;

    inplaceReorder(sortedIndices, surfaces_);
    inplaceReorder(sortedIndices, priorities);
    cellSizeFunctions_.reorder(sortedIndices);

    forAll(priorities, surfI)
    {
        if (priorities[surfI] < defaultPriority_)
        {
            WarningIn("cellSizeControlSurfaces::cellSizeControlSurfaces")
                << "Priority of " << priorities[surfI]
                << " is less than defaultPriority " << defaultPriority_
                << ". All cellSizeFunctions with priorities lower than default "
                << "will be ignored."
                << endl;

            surfaces_.setSize(surfI);
            cellSizeFunctions_.setSize(surfI);

            break;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellSizeControlSurfaces::~cellSizeControlSurfaces()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::cellSizeControlSurfaces::cellSize
(
    const point& pt
) const
{
    scalar size = defaultCellSize_;

    bool refinementTriangulationSwitch = true;

    if (!refinementTriangulationSwitch)
    {
        evalCellSizeFunctions(pt, size);
    }
    else
    {
        Cell_handle ch = T_.locate
        (
            Point(pt.x(), pt.y(), pt.z()),
            oldCellHandle_
        );

        oldCellHandle_ = ch;

        pointFromPoint pA = topoint(ch->vertex(0)->point());
        pointFromPoint pB = topoint(ch->vertex(1)->point());
        pointFromPoint pC = topoint(ch->vertex(2)->point());
        pointFromPoint pD = topoint(ch->vertex(3)->point());

        tetPointRef tet(pA, pB, pC, pD);

        scalarList bary;
        tet.barycentric(pt, bary);

        scalar value = 0;
        forAll(bary, pI)
        {
            value += bary[pI]*ch->vertex(pI)->value();
        }

        size = value;
    }

//if (!anyFunctionFound)
//{
//    // Check if the point in question was actually inside the domain, if
//    // not, then it may be falling back to an inappropriate default size.

//    if (cvMesh_.geometryToConformTo().outside(pt))
//    {
//        pointIndexHit surfHit;
//        label hitSurface;

//        cvMesh_.geometryToConformTo().findSurfaceNearest
//        (
//            pt,
//            sqr(GREAT),
//            surfHit,
//            hitSurface
//        );

//        if (!surfHit.hit())
//        {
//            FatalErrorIn
//            (
//                "Foam::scalar Foam::cellSizeControlSurfaces::cellSize"
//                "("
//                    "const point& pt"
//                ") const"
//            )
//                << "Point " << pt << " did not find a nearest surface point"
//                << nl << exit(FatalError) << endl;
//        }
//    }
//}

    return size;
}


Foam::scalarField Foam::cellSizeControlSurfaces::cellSize
(
    const pointField& pts
) const
{
    scalarField cellSizes(pts.size());

    forAll(pts, i)
    {
        cellSizes[i] = cellSize(pts[i]);
    }

    return cellSizes;
}


void Foam::cellSizeControlSurfaces::setCellSizes
(
    const pointField& pts
)
{
    if (cellSizeFunctions_.size())
    {
        forAll(cellSizeFunctions_, i)
        {
            cellSizeFunction& cSF = cellSizeFunctions_[i];

            cSF.setCellSize(pts);
        }
    }
}


// ************************************************************************* //
