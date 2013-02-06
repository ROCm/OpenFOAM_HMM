/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "controlMeshRefinement.H"
#include "cellSizeAndAlignmentControl.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(controlMeshRefinement, 0);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::controlMeshRefinement::calcFirstDerivative
(
    const Foam::point& a,
    const scalar& cellSizeA,
    const Foam::point& b,
    const scalar& cellSizeB
) const
{
    return (cellSizeA - cellSizeB)/mag(a - b);
}


Foam::scalar Foam::controlMeshRefinement::calcSecondDerivative
(
    const Foam::point& a,
    const scalar& cellSizeA,
    const Foam::point& midPoint,
    const scalar& cellSizeMid,
    const Foam::point& b,
    const scalar& cellSizeB
) const
{
    return (cellSizeA - 2*cellSizeMid + cellSizeB)/magSqr((a - b)/2);
}


bool Foam::controlMeshRefinement::detectEdge
(
    const Foam::point& startPt,
    const Foam::point& endPt,
    pointHit& pointFound,
    const scalar tolSqr,
    const scalar secondDerivTolSqr
) const
{
    Foam::point a(startPt);
    Foam::point b(endPt);

    Foam::point midPoint = (a + b)/2.0;

    label nIterations = 0;

    while (true)
    {
        nIterations++;

        if (magSqr(a - b) < tolSqr)
        {
            pointFound.setPoint(midPoint);
            pointFound.setHit();

            return true;
        }

        // Split into two regions

        const scalar cellSizeA = sizeControls_.cellSize(a);
        const scalar cellSizeB = sizeControls_.cellSize(b);
        const scalar cellSizeMid = sizeControls_.cellSize(midPoint);

        // Region 1
        Foam::point midPoint1 = (a + midPoint)/2.0;
        const scalar cellSizeMid1 = sizeControls_.cellSize(midPoint1);

//        scalar firstDerivative1 =
//            calcFirstDerivative(cellSizeA, cellSizeMid);
        scalar secondDerivative1 =
            calcSecondDerivative
            (
                a,
                cellSizeA,
                midPoint1,
                cellSizeMid1,
                midPoint,
                cellSizeMid
            );

        // Region 2
        Foam::point midPoint2 = (midPoint + b)/2.0;
        const scalar cellSizeMid2 = sizeControls_.cellSize(midPoint2);

//        scalar firstDerivative2 =
//            calcFirstDerivative(f, cellSizeMid, cellSizeB);
        scalar secondDerivative2 =
            calcSecondDerivative
            (
                midPoint,
                cellSizeMid,
                midPoint2,
                cellSizeMid2,
                b,
                cellSizeB
            );

        // Neither region appears to have an inflection
        // To be sure should use higher order derivatives
        if
        (
            magSqr(secondDerivative1) < secondDerivTolSqr
         && magSqr(secondDerivative2) < secondDerivTolSqr
        )
        {
            return false;
        }

        // Pick region with greatest second derivative
        if (magSqr(secondDerivative1) > magSqr(secondDerivative2))
        {
            b = midPoint;
            midPoint = midPoint1;
        }
        else
        {
            a = midPoint;
            midPoint = midPoint2;
        }
    }
}


Foam::pointHit Foam::controlMeshRefinement::findDiscontinuities
(
    const linePointRef& l
) const
{
//    if (true)
//    {
//        Info<< "Searching " << l << " for a discontinuity" << endl;
//    }

    pointHit p(point::max);

    const scalar tolSqr = sqr(1e-6);
    const scalar secondDerivTolSqr = sqr(1e-3);

    detectEdge
    (
        l.start(),
        l.end(),
        p,
        tolSqr,
        secondDerivTolSqr
    );

//    if (p.hit())
//    {
//        Info<< "    " << l << " has a discontinuity at " << p << endl;
//    }
//    else
//    {
//        Info<< "    No discontinuity found" << endl;
//    }

    return p;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::controlMeshRefinement::controlMeshRefinement
(
    cellShapeControl& shapeController
)
:
    shapeController_(shapeController),
    mesh_(shapeController.shapeControlMesh()),
    sizeControls_(shapeController.sizeAndAlignment()),
    geometryToConformTo_(sizeControls_.geometryToConformTo())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::controlMeshRefinement::~controlMeshRefinement()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::controlMeshRefinement::initialMeshPopulation
(
    const autoPtr<backgroundMeshDecomposition>& decomposition
)
{
    autoPtr<boundBox> overallBoundBox;

    // Need to pass in the background mesh decomposition so that can test if
    // a point to insert is on the processor.
    if (Pstream::parRun())
    {
        overallBoundBox.set(new boundBox(decomposition().procBounds()));
    }
    else
    {
        overallBoundBox.set
        (
            new boundBox(geometryToConformTo_.geometry().bounds())
        );
    }

    mesh_.insertBoundingPoints
    (
        overallBoundBox(),
        sizeControls_
    );

    const PtrList<cellSizeAndAlignmentControl>& controlFunctions =
        sizeControls_.controlFunctions();

    forAll(controlFunctions, fI)
    {
        const cellSizeAndAlignmentControl& controlFunction =
            controlFunctions[fI];

        Info<< "Inserting points from " << controlFunction.name()
            << " (" << controlFunction.type() << ")" << endl;

        pointField pts;
        scalarField sizes;
        triadField alignments;

        controlFunction.initialVertices(pts, sizes, alignments);

        List<Vb> vertices(pts.size());

        // Clip the minimum size
        forAll(vertices, vI)
        {
            vertices[vI] = Vb(pts[vI], Vb::vtInternalNearBoundary);
            vertices[vI].targetCellSize() = max
            (
                sizes[vI],
                shapeController_.minimumCellSize()
            );
            vertices[vI].alignment() = alignments[vI];
        }

        pts.clear();
        sizes.clear();
        alignments.clear();

        label nRejected = 0;

        PackedBoolList keepVertex(vertices.size(), true);

        if (Pstream::parRun())
        {
            forAll(vertices, vI)
            {
                const bool onProc = decomposition().positionOnThisProcessor
                (
                    topoint(vertices[vI].point())
                );

                if (!onProc)
                {
                    keepVertex[vI] = false;
                }
            }
        }

        inplaceSubset(keepVertex, vertices);

        const label preInsertedSize = mesh_.number_of_vertices();

        forAll(vertices, vI)
        {
            bool insertPoint = false;

            pointFromPoint pt(topoint(vertices[vI].point()));

            if
            (
                mesh_.dimension() < 3
             || mesh_.is_infinite
                (
                    mesh_.locate(vertices[vI].point())
                )
            )
            {
                insertPoint = true;
            }

            const scalar interpolatedCellSize = shapeController_.cellSize(pt);
            const triad interpolatedAlignment =
                shapeController_.cellAlignment(pt);
            const scalar calculatedCellSize = vertices[vI].targetCellSize();
            const triad calculatedAlignment = vertices[vI].alignment();

            if (debug)
            {
                Info<< "Point = " << pt << nl
                    << "  Size(interp) = " << interpolatedCellSize << nl
                    << "    Size(calc) = " << calculatedCellSize << nl
                    << " Align(interp) = " << interpolatedAlignment << nl
                    << "   Align(calc) = " << calculatedAlignment << nl
                    << endl;
            }

            const scalar sizeDiff =
                mag(interpolatedCellSize - calculatedCellSize);
            const scalar alignmentDiff =
                diff(interpolatedAlignment, calculatedAlignment);

            if (debug)
            {
                Info<< "    size difference = " << sizeDiff
                    << ", alignment difference = " << alignmentDiff << endl;
            }

            // @todo Also need to base it on the alignments
            if
            (
                sizeDiff/interpolatedCellSize > 0.1
             || alignmentDiff > 0.15
            )
            {
                insertPoint = true;
            }

            if (insertPoint)
            {
                mesh_.insert
                (
                    pt,
                    calculatedCellSize,
                    vertices[vI].alignment(),
                    Vb::vtInternalNearBoundary
                );
            }
        }

        //mesh_.rangeInsertWithInfo(vertices.begin(), vertices.end());

        Info<< "    Inserted "
            << returnReduce
               (
                   label(mesh_.number_of_vertices()) - preInsertedSize,
                   sumOp<label>()
               )
            << "/" << returnReduce(vertices.size(), sumOp<label>())
            << endl;
    }
}


Foam::label Foam::controlMeshRefinement::refineMesh
(
    const autoPtr<backgroundMeshDecomposition>& decomposition
)
{
    Info<< "Iterate over cell size mesh edges" << endl;

    DynamicList<Vb> verts(mesh_.number_of_vertices());

    for
    (
        CellSizeDelaunay::Finite_edges_iterator eit =
            mesh_.finite_edges_begin();
        eit != mesh_.finite_edges_end();
        ++eit
    )
    {
        CellSizeDelaunay::Cell_handle c = eit->first;
        CellSizeDelaunay::Vertex_handle vA = c->vertex(eit->second);
        CellSizeDelaunay::Vertex_handle vB = c->vertex(eit->third);

        if
        (
            mesh_.is_infinite(vA)
         || mesh_.is_infinite(vB)
         || (vA->referred() && vB->referred())
         || (vA->referred() && vA->procIndex() > vB->procIndex())
         || (vB->referred() && vB->procIndex() > vA->procIndex())
        )
        {
            continue;
        }

        pointFromPoint ptA(topoint(vA->point()));
        pointFromPoint ptB(topoint(vB->point()));

        linePointRef l(ptA, ptB);

//        Info<< "Edge " << l << endl;

//        Info<< vA->info() << vB->info();

        const Foam::point midPoint = l.centre();
        const scalar length = l.mag();

        scalar evaluatedSize = sizeControls_.cellSize(midPoint);
        scalar interpolatedSize =
            (vA->targetCellSize() + vB->targetCellSize())/2;

        const scalar diff = mag(interpolatedSize - evaluatedSize);

//        Info<< "       Evaluated size = " << evaluatedSize << nl
//            << "    Interpolated size = " << interpolatedSize << nl
//            << "         vA cell size = " << vA->targetCellSize() << nl
//            << "         vB cell size = " << vB->targetCellSize() << endl;

        const pointHit hitPt = findDiscontinuities(l);

        if (hitPt.hit())
        {
            const Foam::point& pt = hitPt.hitPoint();

            if (!geometryToConformTo_.inside(pt))
            {
                continue;
            }

            if (Pstream::parRun())
            {
                if (!decomposition().positionOnThisProcessor(pt))
                {
                    continue;
                }
            }

            verts.append
            (
                Vb
                (
                    toPoint<cellShapeControlMesh::Point>(pt),
                    Vb::vtInternal
                )
            );

            verts.last().targetCellSize() = sizeControls_.cellSize(pt);
            verts.last().alignment() = triad::unset;
        }

//        if (diff/interpolatedSize > 0.1)
//        {
//            // Create a new point
//            // Walk along the edge, placing points where the gradient in cell
//            // size changes
//            if (vA->targetCellSize() <= vB->targetCellSize())
//            {
//                // Walk from vA
//                const label nPoints = length/vA->targetCellSize();
//
//                if (nPoints == 0)
//                {
//                    continue;
//                }
//
//                scalar prevSize = vA->targetCellSize();
//
//                for (label pI = 1; pI < nPoints; ++pI)
//                {
//                    const Foam::point testPt =
//                        topoint(vA->point()) + pI*l.vec()/nPoints;
//
////                    if (geometryToConformTo_.outside(testPt))
////                    {
////                        continue;
////                    }
//
//                    scalar testPtSize = sizeControls_.cellSize(testPt);
//
//                    if (mag(testPtSize - prevSize) < 1e-3*testPtSize)
//                    {
////                        Info<< "Adding " << testPt << " " << testPtSize
////                            << endl;
//                        verts.append
//                        (
//                            Vb
//                            (
//                                toPoint<cellShapeControlMesh::Point>(testPt),
//                                Vb::vtInternal
//                            )
//                        );
//                        verts.last().targetCellSize() = testPtSize;
//                        verts.last().alignment() = triad::unset;
//
//                        break;
//                    }
//
//                    prevSize = testPtSize;
//                }
//            }
//            else
//            {
//                // Walk from vB
//                const label nPoints = length/vB->targetCellSize();
//
//                if (nPoints == 0)
//                {
//                    continue;
//                }
//
//                scalar prevSize = vA->targetCellSize();
//
//                for (label pI = 1; pI < nPoints; ++pI)
//                {
//                    const Foam::point testPt =
//                        topoint(vB->point()) - pI*l.vec()/nPoints;
//
////                    if (geometryToConformTo_.outside(testPt))
////                    {
////                        continue;
////                    }
//
//                    scalar testPtSize = sizeControls_.cellSize(testPt);
//
//                    if (mag(testPtSize - prevSize) < 1e-3*testPtSize)
//                    //if (testPtSize > prevSize)// || testPtSize < prevSize)
//                    {
////                        Info<< "Adding " << testPt << " " << testPtSize
////                            << endl;
//                        verts.append
//                        (
//                            Vb
//                            (
//                                toPoint<cellShapeControlMesh::Point>(testPt),
//                                Vb::vtInternal
//                            )
//                        );
//                        verts.last().targetCellSize() = testPtSize;
//                        verts.last().alignment() = triad::unset;
//
//                        break;
//                    }
//
//                    prevSize = testPtSize;
//                }
//            }
//        }
    }

//    const pointField cellCentres(mesh_.cellCentres());
//
//    Info<< "    Created cell centres" << endl;
//
//    const PtrList<cellSizeAndAlignmentControl>& controlFunctions =
//        sizeControls_.controlFunctions();
//
//    DynamicList<Vb> verts(mesh_.number_of_vertices());
//
//    forAll(cellCentres, cellI)
//    {
//        Foam::point pt = cellCentres[cellI];
//
//        if (!geometryToConformTo_.inside(pt))
//        {
//            continue;
//        }
//
//        scalarList bary;
//        cellShapeControlMesh::Cell_handle ch;
//
//        mesh_.barycentricCoords(pt, bary, ch);
//
//        if (mesh_.is_infinite(ch))
//        {
//            continue;
//        }
//
//        scalar interpolatedSize = 0;
//        forAll(bary, pI)
//        {
//            interpolatedSize += bary[pI]*ch->vertex(pI)->targetCellSize();
//        }
//
//        label lastPriority = labelMin;
//        scalar lastCellSize = GREAT;
//        forAll(controlFunctions, fI)
//        {
//            const cellSizeAndAlignmentControl& controlFunction =
//                controlFunctions[fI];
//
//            if (controlFunction.priority() < lastPriority)
//            {
//                continue;
//            }
//
//            if (isA<searchableSurfaceControl>(controlFunction))
//            {
//                const cellSizeFunction& sizeFunction =
//                    dynamicCast<const searchableSurfaceControl>
//                    (
//                        controlFunction
//                    ).sizeFunction();
//
//                scalar cellSize = 0;
//                if (sizeFunction.cellSize(pt, cellSize))
//                {
//                    if (controlFunction.priority() == lastPriority)
//                    {
//                        if (cellSize < lastCellSize)
//                        {
//                            lastCellSize = cellSize;
//                        }
//                    }
//                    else
//                    {
//                        lastCellSize = cellSize;
//                    }
//
//                    lastPriority = controlFunction.priority();
//                }
//            }
//        }
//
//        const scalar diff = mag(interpolatedSize - lastCellSize);
//
//        if (diff/interpolatedSize > 0.1)
//        {
//            bool pointFound = false;
//
//            // Travel along lines to each tet point
//            for (label vI = 0; vI < 4; ++vI)
//            {
//                if (pointFound)
//                {
//                    break;
//                }
//
//                pointFromPoint tetVertex = topoint(ch->vertex(vI)->point());
//                scalar tetVertexSize = sizeControls_.cellSize(tetVertex);
//
//                if (tetVertexSize < interpolatedSize)
//                {
//                    linePointRef l(tetVertex, pt);
//
//                    label nTestPoints = l.mag()/tetVertexSize;
//
//                    for (label i = 1; i < nTestPoints; ++i)
//                    {
//                        const Foam::point testPt
//                        (
//                            tetVertex
//                          + l.vec()*i/nTestPoints
//                        );
//
//                        scalar testPtSize = sizeControls_.cellSize(testPt);
//
//                        if (testPtSize > tetVertexSize)
//                        {
//                            pt = testPt;
//                            lastCellSize = testPtSize;
//                            pointFound = true;
//                            break;
//                        }
//                    }
//                }
//                else
//                {
//                    linePointRef l(pt, tetVertex);
//
//                    label nTestPoints = l.mag()/tetVertexSize;
//
//                    for (label i = 1; i < nTestPoints; ++i)
//                    {
//                        const Foam::point testPt
//                        (
//                            tetVertex
//                          + l.vec()*i/nTestPoints
//                        );
//
//                        scalar testPtSize = sizeControls_.cellSize(testPt);
//
//                        if (testPtSize < interpolatedSize)
//                        {
//                            pt = testPt;
//                            lastCellSize = testPtSize;
//                            pointFound = true;
//                            break;
//                        }
//                    }
//                }
//            }
//        }
//
//        // Decide whether to insert or not
//        if (lastCellSize < GREAT)
//        {
//            if (!geometryToConformTo_.inside(pt))
//            {
//                continue;
//            }
//
//            if (Pstream::parRun())
//            {
//                if (!decomposition().positionOnThisProcessor(pt))
//                {
//                    continue;
//                }
//            }
//
//            verts.append
//            (
//                Vb
//                (
//                    toPoint<cellShapeControlMesh::Point>(pt),
//                    Vb::vtInternal
//                )
//            );
//            verts.last().targetCellSize() = lastCellSize;
//            verts.last().alignment() = triad::unset;
//        }
//    }

    mesh_.insertPoints(verts);

    return verts.size();
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
