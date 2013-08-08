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

#include "cellShapeControlMesh.H"
#include "cellSizeAndAlignmentControls.H"
#include "pointIOField.H"
#include "scalarIOField.H"
#include "triadIOField.H"
#include "tetrahedron.H"
#include "plane.H"
#include "transform.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(cellShapeControlMesh, 0);

word cellShapeControlMesh::meshSubDir = "cellShapeControlMesh";
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

//Foam::tensor Foam::cellShapeControlMesh::requiredAlignment
//(
//    const Foam::point& pt,
//    const searchableSurfaces& allGeometry,
//    const conformationSurfaces& geometryToConformTo
//) const
//{
//    pointIndexHit surfHit;
//    label hitSurface;
//
//    geometryToConformTo.findSurfaceNearest
//    (
//        pt,
//        sqr(GREAT),
//        surfHit,
//        hitSurface
//    );
//
//    if (!surfHit.hit())
//    {
//        FatalErrorIn
//        (
//            "Foam::tensor Foam::conformalVoronoiMesh::requiredAlignment"
//        )   << "findSurfaceNearest did not find a hit across the surfaces."
//            << exit(FatalError) << endl;
//    }
//
//    // Primary alignment
//
//    vectorField norm(1);
//
//    allGeometry[hitSurface].getNormal
//    (
//        List<pointIndexHit>(1, surfHit),
//        norm
//    );
//
//    const vector np = norm[0];
//
//    // Generate equally spaced 'spokes' in a circle normal to the
//    // direction from the vertex to the closest point on the surface
//    // and look for a secondary intersection.
//
//    const vector d = surfHit.hitPoint() - pt;
//
//    const tensor Rp = rotationTensor(vector(0,0,1), np);
//
//    const label s = 36;//foamyHexMeshControls().alignmentSearchSpokes();
//
//    scalar closestSpokeHitDistance = GREAT;
//
//    pointIndexHit closestSpokeHit;
//
//    label closestSpokeSurface = -1;
//
//    const scalar spanMag = geometryToConformTo.globalBounds().mag();
//
//    for (label i = 0; i < s; i++)
//    {
//        vector spoke
//        (
//            Foam::cos(i*constant::mathematical::twoPi/s),
//            Foam::sin(i*constant::mathematical::twoPi/s),
//            0
//        );
//
//        spoke *= spanMag;
//
//        spoke = Rp & spoke;
//
//        pointIndexHit spokeHit;
//
//        label spokeSurface = -1;
//
//        // internal spoke
//
//        geometryToConformTo.findSurfaceNearestIntersection
//        (
//            pt,
//            pt + spoke,
//            spokeHit,
//            spokeSurface
//        );
//
//        if (spokeHit.hit())
//        {
//            scalar spokeHitDistance = mag
//            (
//                spokeHit.hitPoint() - pt
//            );
//
//            if (spokeHitDistance < closestSpokeHitDistance)
//            {
//                closestSpokeHit = spokeHit;
//                closestSpokeSurface = spokeSurface;
//                closestSpokeHitDistance = spokeHitDistance;
//            }
//        }
//
//        //external spoke
//
//        Foam::point mirrorPt = pt + 2*d;
//
//        geometryToConformTo.findSurfaceNearestIntersection
//        (
//            mirrorPt,
//            mirrorPt + spoke,
//            spokeHit,
//            spokeSurface
//        );
//
//        if (spokeHit.hit())
//        {
//            scalar spokeHitDistance = mag
//            (
//                spokeHit.hitPoint() - mirrorPt
//            );
//
//            if (spokeHitDistance < closestSpokeHitDistance)
//            {
//                closestSpokeHit = spokeHit;
//                closestSpokeSurface = spokeSurface;
//                closestSpokeHitDistance = spokeHitDistance;
//            }
//        }
//    }
//
//    if (closestSpokeSurface == -1)
//    {
////        WarningIn
////        (
////            "conformalVoronoiMesh::requiredAlignment"
////            "("
////                "const Foam::point& pt"
////            ") const"
////        )   << "No secondary surface hit found in spoke search "
////            << "using " << s
////            << " spokes, try increasing alignmentSearchSpokes."
////            << endl;
//
//        return I;
//    }
//
//    // Auxiliary alignment generated by spoke intersection normal.
//
//    allGeometry[closestSpokeSurface].getNormal
//    (
//        List<pointIndexHit>(1, closestSpokeHit),
//        norm
//    );
//
//    const vector& na = norm[0];
//
//    // Secondary alignment
//    vector ns = np ^ na;
//
//    if (mag(ns) < SMALL)
//    {
//        FatalErrorIn("conformalVoronoiMesh::requiredAlignment")
//            << "Parallel normals detected in spoke search." << nl
//            << "point: " << pt << nl
//            << "closest surface point: " << surfHit.hitPoint() << nl
//            << "closest spoke hit: " << closestSpokeHit.hitPoint() << nl
//            << "np: " << surfHit.hitPoint() + np << nl
//            << "ns: " << closestSpokeHit.hitPoint() + na << nl
//            << exit(FatalError);
//    }
//
//    ns /= mag(ns);
//
//    tensor Rs = rotationTensor((Rp & vector(0,1,0)), ns);
//
//    return (Rs & Rp);
//}


Foam::label Foam::cellShapeControlMesh::removePoints()
{
    label nRemoved = 0;
    for
    (
        CellSizeDelaunay::Finite_vertices_iterator vit =
            finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        std::list<Vertex_handle> verts;
        adjacent_vertices(vit, std::back_inserter(verts));

        bool removePt = true;
        for
        (
            std::list<Vertex_handle>::iterator aVit = verts.begin();
            aVit != verts.end();
            ++aVit
        )
        {
            Vertex_handle avh = *aVit;

            scalar diff =
                mag(avh->targetCellSize() - vit->targetCellSize())
               /max(vit->targetCellSize(), 1e-6);

            if (diff > 0.05)
            {
                removePt = false;
            }
        }

        if (removePt)
        {
            remove(vit);
            nRemoved++;
        }
    }

    return nRemoved;
}


Foam::tmp<Foam::pointField> Foam::cellShapeControlMesh::cellCentres() const
{
    tmp<pointField> tcellCentres(new pointField(number_of_finite_cells()));
    pointField& cellCentres = tcellCentres();

    label count = 0;
    for
    (
        CellSizeDelaunay::Finite_cells_iterator c = finite_cells_begin();
        c != finite_cells_end();
        ++c
    )
    {
        if (c->hasFarPoint())
        {
            continue;
        }

        scalarList bary;
        cellShapeControlMesh::Cell_handle ch;

        const Foam::point centre = topoint
        (
            CGAL::centroid<baseK>
            (
                c->vertex(0)->point(),
                c->vertex(1)->point(),
                c->vertex(2)->point(),
                c->vertex(3)->point()
            )
        );

        cellCentres[count++] = centre;
    }

    cellCentres.resize(count);

    return tcellCentres;
}


void Foam::cellShapeControlMesh::writeTriangulation()
{
    OFstream str
    (
        "refinementTriangulation_"
      + name(Pstream::myProcNo())
      + ".obj"
    );

    label count = 0;

    Info<< "Write refinementTriangulation" << endl;

    for
    (
        CellSizeDelaunay::Finite_edges_iterator e = finite_edges_begin();
        e != finite_edges_end();
        ++e
    )
    {
        Cell_handle c = e->first;
        Vertex_handle vA = c->vertex(e->second);
        Vertex_handle vB = c->vertex(e->third);

        // Don't write far edges
        if (vA->farPoint() || vB->farPoint())
        {
            continue;
        }

        // Don't write unowned edges
        if (vA->referred() && vB->referred())
        {
            continue;
        }

        pointFromPoint p1 = topoint(vA->point());
        pointFromPoint p2 = topoint(vB->point());

        meshTools::writeOBJ(str, p1, p2, count);
    }

    if (is_valid())
    {
        Info<< "    Triangulation is valid" << endl;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::triangulatedMesh::writeRefinementTriangulation()"
        )   << "Triangulation is not valid"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellShapeControlMesh::cellShapeControlMesh(const Time& runTime)
:
    DistributedDelaunayMesh<CellSizeDelaunay>
    (
        runTime,
        meshSubDir
    ),
    runTime_(runTime),
    defaultCellSize_(0.0)
{
    if (this->vertexCount())
    {
        fvMesh mesh
        (
            IOobject
            (
                meshSubDir,
                runTime.timeName(),
                runTime,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        );

        if (mesh.nPoints() == this->vertexCount())
        {
            pointScalarField sizes
            (
                IOobject
                (
                    "sizes",
                    runTime.timeName(),
                    meshSubDir,
                    runTime,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                pointMesh::New(mesh)
            );

            triadIOField alignments
            (
                IOobject
                (
                    "alignments",
                    mesh.time().timeName(),
                    meshSubDir,
                    mesh.time(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                )
            );

            if
            (
                sizes.size() == this->vertexCount()
             && alignments.size() == this->vertexCount()
            )
            {
                for
                (
                    Finite_vertices_iterator vit = finite_vertices_begin();
                    vit != finite_vertices_end();
                    ++vit
                )
                {
                    vit->targetCellSize() = sizes[vit->index()];
                    vit->alignment() = alignments[vit->index()];
                }
            }
            else
            {
                FatalErrorIn
                (
                    "Foam::cellShapeControlMesh::cellShapeControlMesh"
                    "(const Time&)"
                )   << "Cell size point field is not the same size as the "
                    << "mesh."
                    << abort(FatalError);
            }
        }
    }
}


//Foam::triangulatedMesh::triangulatedMesh
//(
//    const Time& runTime,
//    const fileName& pointsFile,
//    const fileName& sizesFile,
//    const fileName& alignmentsFile,
//    const scalar& defaultCellSize
//)
//:
//    defaultCellSize_(defaultCellSize)
//{
//    Info<< "    Reading points from file     : " << pointsFile << endl;
//
//    pointIOField points
//    (
//        IOobject
//        (
//            pointsFile,
//            runTime.constant(),
//            runTime,
//            IOobject::MUST_READ,
//            IOobject::NO_WRITE
//        )
//    );
//
//    Info<< "    Reading sizes from file      : " << sizesFile << endl;
//
//    scalarIOField sizes
//    (
//        IOobject
//        (
//            sizesFile,
//            runTime.constant(),
//            runTime,
//            IOobject::MUST_READ,
//            IOobject::NO_WRITE
//        )
//    );
//
//    Info<< "    Reading alignments from file : " << alignmentsFile << endl;
//
//    tensorIOField alignments
//    (
//        IOobject
//        (
//            alignmentsFile,
//            runTime.constant(),
//            runTime,
//            IOobject::MUST_READ,
//            IOobject::NO_WRITE
//        )
//    );
//
//    Info<< "    Number of points : " << points.size() << endl;
//    Info<< "    Minimum size     : " << min(sizes) << endl;
//    Info<< "    Average size     : " << average(sizes) << endl;
//    Info<< "    Maximum size     : " << max(sizes) << endl;
//
//    forAll(points, pI)
//    {
//        size_t nVert = number_of_vertices();
//
//        Vertex_handle v = insert
//        (
//            Point(points[pI].x(), points[pI].y(), points[pI].z())
//        );
//
//        if (number_of_vertices() == nVert)
//        {
//            Info<< "    Failed to insert point : " << points[pI] << endl;
//        }
//
//        v->targetCellSize() = sizes[pI];
//
//        const tensor& alignment = alignments[pI];
//
//
//
//        v->alignment() = alignment;
//    }
//
////    scalar factor = 1.0;
////    label maxIteration = 1;
////
////    for (label iteration = 0; iteration < maxIteration; ++iteration)
////    {
////        Info<< "Iteration : " << iteration << endl;
////
////        label nRefined = refineTriangulation(factor);
////
////        Info<< "    Number of cells refined in refinement iteration : "
////            << nRefined << nl << endl;
////
////        if (nRefined <= 0 && iteration != 0)
////        {
////            break;
////        }
////
////        factor *= 1.5;
////    }
//
//    //writeRefinementTriangulation();
//}


//Foam::triangulatedMesh::triangulatedMesh
//(
//    const Time& runTime,
//    const DynamicList<Foam::point>& points,
//    const DynamicList<scalar>& sizes,
//    const DynamicList<tensor>& alignments,
//    const scalar& defaultCellSize
//)
//:
//    defaultCellSize_(defaultCellSize)
//{
//    forAll(points, pI)
//    {
//        size_t nVert = number_of_vertices();
//
//        Vertex_handle v = insert
//        (
//            Point(points[pI].x(), points[pI].y(), points[pI].z())
//        );
//
//        if (number_of_vertices() == nVert)
//        {
//            Info<< "Failed to insert point : " << points[pI] << endl;
//        }
//
//        v->targetCellSize() = sizes[pI];
//
//        v->alignment() = alignments[pI];
//    }
//
//    //writeRefinementTriangulation();
//
//    Info<< nl << "Refinement triangulation information: " << endl;
//    Info<< "    Number of vertices: " << label(number_of_vertices()) << endl;
//    Info<< "    Number of cells   : "
//        << label(number_of_finite_cells()) << endl;
//    Info<< "    Number of faces   : "
//        << label(number_of_finite_facets()) << endl;
//    Info<< "    Number of edges   : "
//        << label(number_of_finite_edges()) << endl;
//    Info<< "    Dimensionality    : " << label(dimension()) << nl << endl;
//}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellShapeControlMesh::~cellShapeControlMesh()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::cellShapeControlMesh::barycentricCoords
(
    const Foam::point& pt,
    scalarList& bary,
    Cell_handle& ch
) const
{
    // Use the previous cell handle as a hint on where to start searching
    // Giving a hint causes strange errors...
    ch = locate(toPoint(pt));

    if (dimension() > 2 && !is_infinite(ch))
    {
        oldCellHandle_ = ch;

        tetPointRef tet
        (
            topoint(ch->vertex(0)->point()),
            topoint(ch->vertex(1)->point()),
            topoint(ch->vertex(2)->point()),
            topoint(ch->vertex(3)->point())
        );

        tet.barycentric(pt, bary);
    }
}


Foam::boundBox Foam::cellShapeControlMesh::bounds() const
{
    DynamicList<Foam::point> pts(number_of_vertices());

    for
    (
        Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->real())
        {
            pts.append(topoint(vit->point()));
        }
    }

    boundBox bb(pts);

    return bb;
}


void Foam::cellShapeControlMesh::distribute
(
    const backgroundMeshDecomposition& decomposition
)
{
    DynamicList<Foam::point> points(number_of_vertices());
    DynamicList<scalar> sizes(number_of_vertices());
    DynamicList<tensor> alignments(number_of_vertices());

    DynamicList<Vb> farPts(8);

    for
    (
        Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->real())
        {
            points.append(topoint(vit->point()));
            sizes.append(vit->targetCellSize());
            alignments.append(vit->alignment());
        }
        else if (vit->farPoint())
        {
            farPts.append
            (
                Vb
                (
                    vit->point(),
                    -1,
                    Vb::vtFar,
                    Pstream::myProcNo()
                )
            );

            farPts.last().targetCellSize() = vit->targetCellSize();
            farPts.last().alignment() = vit->alignment();
        }
    }

    autoPtr<mapDistribute> mapDist =
        DistributedDelaunayMesh<CellSizeDelaunay>::distribute
        (
            decomposition,
            points
        );

    mapDist().distribute(sizes);
    mapDist().distribute(alignments);

    // Reset the entire tessellation
    DelaunayMesh<CellSizeDelaunay>::reset();


    // Internal points have to be inserted first
    DynamicList<Vb> verticesToInsert(points.size());


    forAll(farPts, ptI)
    {
        verticesToInsert.append(farPts[ptI]);
    }


    forAll(points, pI)
    {
        verticesToInsert.append
        (
            Vb
            (
                toPoint(points[pI]),
                -1,
                Vb::vtInternal,
                Pstream::myProcNo()
            )
        );

        verticesToInsert.last().targetCellSize() = sizes[pI];
        verticesToInsert.last().alignment() = alignments[pI];
    }

    Info<< nl << "    Inserting distributed background tessellation..." << endl;

    this->rangeInsertWithInfo
    (
        verticesToInsert.begin(),
        verticesToInsert.end(),
        true
    );

    sync(decomposition.procBounds());

    Info<< "    Total number of vertices after redistribution "
        << returnReduce(label(number_of_vertices()), sumOp<label>()) << endl;
}


Foam::tensorField Foam::cellShapeControlMesh::dumpAlignments() const
{
    tensorField alignmentsTmp(number_of_vertices(), tensor::zero);

    label count = 0;
    for
    (
        Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        alignmentsTmp[count++] = vit->alignment();
    }

    return alignmentsTmp;
}


void Foam::cellShapeControlMesh::insertBoundingPoints
(
    const boundBox& bb,
    const cellSizeAndAlignmentControls& sizeControls
)
{
    // Loop over bound box points and get cell size and alignment
    const pointField bbPoints(bb.points());

    forAll(bbPoints, pI)
    {
        const Foam::point& pt = bbPoints[pI];

        // Cell size here will return default cell size
        const scalar cellSize = sizeControls.cellSize(pt);

        if (debug)
        {
            Info<< "Insert Bounding Point: " << pt << " " << cellSize << endl;
        }

        // Get the cell size of the nearest surface.
//        geometryToConformTo_.findSurfaceNearest
//        (
//            pt,
//            GREAT,
//            surfHit,
//            hitSurface
//        );

        const tensor alignment = tensor::I;

        insert
        (
            pt,
            cellSize,
            alignment,
            Vb::vtInternalNearBoundary
        );
    }
}


void Foam::cellShapeControlMesh::write() const
{
    Info<< "Writing " << meshSubDir << endl;

    // Reindex the cells
    label cellCount = 0;
    for
    (
        Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        if (!cit->hasFarPoint() && !is_infinite(cit))
        {
            cit->cellIndex() = cellCount++;
        }
    }

    DelaunayMesh<CellSizeDelaunay>::labelTolabelPairHashTable vertexMap;
    labelList cellMap;

    autoPtr<polyMesh> meshPtr = DelaunayMesh<CellSizeDelaunay>::createMesh
    (
        meshSubDir,
        vertexMap,
        cellMap
    );
    const polyMesh& mesh = meshPtr();

    pointScalarField sizes
    (
        IOobject
        (
            "sizes",
            mesh.time().timeName(),
            meshSubDir,
            mesh.time(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(mesh),
        scalar(0)
    );

    triadIOField alignments
    (
        IOobject
        (
            "alignments",
            mesh.time().timeName(),
            meshSubDir,
            mesh.time(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sizes.size()
    );

    // Write alignments
//    OFstream str(runTime_.path()/"alignments.obj");

    for
    (
        Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (!vit->farPoint())
        {
            // Populate sizes
            sizes[vertexMap[labelPair(vit->index(), vit->procIndex())]] =
                vit->targetCellSize();

            alignments[vertexMap[labelPair(vit->index(), vit->procIndex())]] =
                vit->alignment();

//            // Write alignments
//            const tensor& alignment = vit->alignment();
//            pointFromPoint pt = topoint(vit->point());
//
//            if
//            (
//                alignment.x() == triad::unset[0]
//             || alignment.y() == triad::unset[0]
//             || alignment.z() == triad::unset[0]
//            )
//            {
//                Info<< "Bad alignment = " << vit->info();
//
//                vit->alignment() = tensor::I;
//
//                Info<< "New alignment = " << vit->info();
//
//                continue;
//            }
//
//            meshTools::writeOBJ(str, pt, alignment.x() + pt);
//            meshTools::writeOBJ(str, pt, alignment.y() + pt);
//            meshTools::writeOBJ(str, pt, alignment.z() + pt);
        }
    }

    mesh.write();
    sizes.write();
    alignments.write();
}


// ************************************************************************* //
