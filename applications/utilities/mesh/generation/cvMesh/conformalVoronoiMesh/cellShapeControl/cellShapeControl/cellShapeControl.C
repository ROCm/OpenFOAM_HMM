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

#include "cellShapeControl.H"
#include "pointField.H"
#include "scalarField.H"
#include "cellSizeAndAlignmentControl.H"
#include "searchableSurfaceControl.H"
#include "cellSizeFunction.H"
#include "triad.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template <class Triangulation, class Type>
Foam::tmp<Foam::Field<Type> > Foam::cellShapeControl::filterFarPoints
(
    const Triangulation& mesh,
    const Field<Type>& field
)
{
    tmp<Field<Type> > tNewField(new Field<Type>(field.size()));
    Field<Type>& newField = tNewField();

    label added = 0;
    label count = 0;

    for
    (
        typename Triangulation::Finite_vertices_iterator vit =
            mesh.finite_vertices_begin();
        vit != mesh.finite_vertices_end();
        ++vit
    )
    {
        if (vit->real())
        {
            newField[added++] = field[count];
        }

        count++;
    }

    newField.resize(added);

    return tNewField;
}


template <class Triangulation>
Foam::autoPtr<Foam::mapDistribute> Foam::cellShapeControl::buildReferredMap
(
    const Triangulation& mesh,
    labelList& indices
)
{
    globalIndex globalIndexing(mesh.vertexCount());

    DynamicList<label> dynIndices(mesh.vertexCount()/10);

    for
    (
        typename Triangulation::Finite_vertices_iterator vit =
            mesh.finite_vertices_begin();
        vit != mesh.finite_vertices_end();
        ++vit
    )
    {
        if (vit->referred())
        {
            dynIndices.append
            (
                globalIndexing.toGlobal(vit->procIndex(), vit->index())
            );
        }
    }

    indices.transfer(dynIndices);

    List<Map<label> > compactMap;
    return autoPtr<mapDistribute>
    (
        new mapDistribute
        (
            globalIndexing,
            indices,
            compactMap
        )
    );
}


template <class Triangulation>
Foam::autoPtr<Foam::mapDistribute> Foam::cellShapeControl::buildMap
(
    const Triangulation& mesh,
    labelListList& pointPoints
)
{
    pointPoints.setSize(mesh.vertexCount());

    globalIndex globalIndexing(mesh.vertexCount());

    for
    (
        typename Triangulation::Finite_vertices_iterator vit =
            mesh.finite_vertices_begin();
        vit != mesh.finite_vertices_end();
        ++vit
    )
    {
        if (!vit->real())
        {
            continue;
        }

        std::list<typename Triangulation::Vertex_handle> adjVerts;
        mesh.finite_adjacent_vertices(vit, std::back_inserter(adjVerts));

        DynamicList<label> indices(adjVerts.size());

        for
        (
            typename std::list<typename Triangulation::Vertex_handle>::
                const_iterator adjVertI = adjVerts.begin();
            adjVertI != adjVerts.end();
            ++adjVertI
        )
        {
            typename Triangulation::Vertex_handle vh = *adjVertI;

            if (!vh->farPoint())
            {
                indices.append
                (
                    globalIndexing.toGlobal(vh->procIndex(), vh->index())
                );
            }
        }

        pointPoints[vit->index()].transfer(indices);
    }

    List<Map<label> > compactMap;
    return autoPtr<mapDistribute>
    (
        new mapDistribute
        (
            globalIndexing,
            pointPoints,
            compactMap
        )
    );
}


template <class Triangulation>
Foam::tmp<Foam::triadField> Foam::cellShapeControl::buildAlignmentField
(
    const Triangulation& mesh
)
{
    tmp<triadField> tAlignments
    (
        new triadField(mesh.vertexCount(), triad::unset)
    );
    triadField& alignments = tAlignments();

    for
    (
        typename Triangulation::Finite_vertices_iterator vit =
            mesh.finite_vertices_begin();
        vit != mesh.finite_vertices_end();
        ++vit
    )
    {
        if (!vit->real())
        {
            continue;
        }

        alignments[vit->index()] = triad
        (
            vit->alignment().x(),
            vit->alignment().y(),
            vit->alignment().z()
        );
    }

    return tAlignments;
}


template <class Triangulation>
Foam::tmp<Foam::pointField> Foam::cellShapeControl::buildPointField
(
    const Triangulation& mesh
)
{
    tmp<pointField> tPoints
    (
        new pointField(mesh.vertexCount(), point(GREAT, GREAT, GREAT))
    );
    pointField& points = tPoints();

    for
    (
        typename Triangulation::Finite_vertices_iterator vit =
            mesh.finite_vertices_begin();
        vit != mesh.finite_vertices_end();
        ++vit
    )
    {
        if (!vit->real())
        {
            continue;
        }

        points[vit->index()] = topoint(vit->point());
    }

    return tPoints;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellShapeControl::cellShapeControl
(
    const Time& runTime,
    const dictionary& motionDict,
    const searchableSurfaces& allGeometry,
    const conformationSurfaces& geometryToConformTo
)
:
    dictionary(motionDict),
    runTime_(runTime),
    allGeometry_(allGeometry),
    geometryToConformTo_(geometryToConformTo),
    defaultCellSize_(readScalar(lookup("defaultCellSize"))),
    shapeControlMesh_(runTime),
    aspectRatio_(motionDict),
    sizeAndAlignment_
    (
        runTime,
        motionDict.subDict("shapeControlFunctions"),
        geometryToConformTo
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellShapeControl::~cellShapeControl()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalarField Foam::cellShapeControl::cellSize
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


Foam::scalar Foam::cellShapeControl::cellSize(const point& pt) const
{
    scalarList bary;
    cellShapeControlMesh::Cell_handle ch;

    shapeControlMesh_.barycentricCoords(pt, bary, ch);

    scalar size = 0;

    label nFarPoints = 0;
    for (label pI = 0; pI < 4; ++pI)
    {
        if (ch->vertex(pI)->farPoint())
        {
            nFarPoints++;
        }
    }

    if (shapeControlMesh_.is_infinite(ch))
    {
//        if (nFarPoints != 0)
//        {
//            for (label pI = 0; pI < 4; ++pI)
//            {
//                if (!ch->vertex(pI)->farPoint())
//                {
//                    size = ch->vertex(pI)->targetCellSize();
//                    return size;
//                }
//            }
//        }

        // Look up nearest point
        cellShapeControlMesh::Vertex_handle nearV =
            shapeControlMesh_.nearest_vertex
            (
                toPoint<cellShapeControlMesh::Point>(pt)
            );

        size = nearV->targetCellSize();
    }
    else
    {
        if (nFarPoints != 0)
        {
            for (label pI = 0; pI < 4; ++pI)
            {
                if (!ch->vertex(pI)->farPoint())
                {
                    size = ch->vertex(pI)->targetCellSize();
                    return size;
                }
            }
        }
        else
        {
            forAll(bary, pI)
            {
                size += bary[pI]*ch->vertex(pI)->targetCellSize();
            }
        }
    }

    return size;
}


//- Return the cell alignment at the given location
Foam::tensor Foam::cellShapeControl::cellAlignment(const point& pt) const
{
    scalarList bary;
    cellShapeControlMesh::Cell_handle ch;

    shapeControlMesh_.barycentricCoords(pt, bary, ch);

    tensor alignment = tensor::zero;

    label nFarPoints = 0;
    for (label pI = 0; pI < 4; ++pI)
    {
        if (ch->vertex(pI)->farPoint())
        {
            nFarPoints++;
        }
    }

    if (shapeControlMesh_.is_infinite(ch) || nFarPoints == 4)
    {
        Pout<< "At Infinite vertex" << endl;

        if (nFarPoints != 0)
        {
            for (label pI = 0; pI < 4; ++pI)
            {
                if (!ch->vertex(pI)->farPoint())
                {
                    alignment = ch->vertex(pI)->alignment();
                    return alignment;
                }
            }
        }

//        cellShapeControlMesh::Vertex_handle nearV =
//            shapeControlMesh_.nearest_vertex
//            (
//                toPoint<cellShapeControlMesh::Point>(pt)
//            );
//
//        alignment = nearV->alignment();
    }
    else
    {
//        forAll(bary, pI)
//        {
//            alignment += bary[pI]*ch->vertex(pI)->alignment();
//        }

        cellShapeControlMesh::Vertex_handle nearV =
            shapeControlMesh_.nearest_vertex_in_cell
            (
                toPoint<cellShapeControlMesh::Point>(pt),
                ch
            );

        alignment = nearV->alignment();
    }

    return alignment;
}


void Foam::cellShapeControl::cellSizeAndAlignment
(
    const point& pt,
    scalar& size,
    tensor& alignment
) const
{
    scalarList bary;
    cellShapeControlMesh::Cell_handle ch;

    shapeControlMesh_.barycentricCoords(pt, bary, ch);

    alignment = tensor::zero;
    size = 0;

    label nFarPoints = 0;
    for (label pI = 0; pI < 4; ++pI)
    {
        if (ch->vertex(pI)->farPoint())
        {
            nFarPoints++;
        }
    }

    if (shapeControlMesh_.is_infinite(ch))
    {
        if (nFarPoints != 0)
        {
            for (label pI = 0; pI < 4; ++pI)
            {
                if (!ch->vertex(pI)->farPoint())
                {
                    size = ch->vertex(pI)->targetCellSize();
                    alignment = ch->vertex(pI)->alignment();
                    return;
                }
            }
        }

//        cellShapeControlMesh::Vertex_handle nearV =
//            shapeControlMesh_.nearest_vertex
//            (
//                toPoint<cellShapeControlMesh::Point>(pt)
//            );
//
//        size = nearV->targetCellSize();
//        alignment = nearV->alignment();
    }
    else
    {
        if (nFarPoints != 0)
        {
            for (label pI = 0; pI < 4; ++pI)
            {
                if (!ch->vertex(pI)->farPoint())
                {
                    size = ch->vertex(pI)->targetCellSize();
                    alignment = ch->vertex(pI)->alignment();
                    return;
                }
            }
        }
        else
        {
            triad tri;

            forAll(bary, pI)
            {
                size += bary[pI]*ch->vertex(pI)->targetCellSize();

                triad triTmp2
                (
                    ch->vertex(pI)->alignment().x(),
                    ch->vertex(pI)->alignment().y(),
                    ch->vertex(pI)->alignment().z()
                );

                tri += triTmp2*bary[pI];
            }

            tri.normalize();
            tri.orthogonalize();
            tri = tri.sortxyz();

            alignment = tensor
            (
                tri.x(), tri.y(), tri.z()
            );
        }
    }
}


void Foam::cellShapeControl::initialMeshPopulation
(
    const autoPtr<backgroundMeshDecomposition>& decomposition
)
{
    // Need to pass in the background mesh decomposition so that can test if
    // a point to insert is on the processor.
    if (Pstream::parRun())
    {
        shapeControlMesh_.insertBoundingPoints(decomposition().procBounds());
    }
    else
    {
        shapeControlMesh_.insertBoundingPoints(allGeometry_.bounds());
    }

    const PtrList<cellSizeAndAlignmentControl>& controlFunctions =
        sizeAndAlignment_.controlFunctions();

    forAll(controlFunctions, fI)
    {
        const cellSizeAndAlignmentControl& controlFunction =
            controlFunctions[fI];

        Info<< "Inserting points from " << controlFunction.name()
            << " (" << controlFunction.type() << ")" << endl;

        pointField pts(1);
        scalarField sizes(1);
        Field<triad> alignments(1);

        controlFunction.initialVertices(pts, sizes, alignments);

        label nRejected = 0;

        forAll(pts, pI)
        {
            if (Pstream::parRun())
            {
                if (!decomposition().positionOnThisProcessor(pts[pI]))
                {
                    nRejected++;
                    continue;
                }
            }

            shapeControlMesh_.insert
            (
                pts[pI],
                sizes[pI],
                alignments[pI]
            );
        }

        Info<< "    Inserted " << (pts.size() - nRejected) << "/" << pts.size()
            << endl;
    }
}


Foam::label Foam::cellShapeControl::refineMesh
(
    const autoPtr<backgroundMeshDecomposition>& decomposition
)
{
    const pointField cellCentres = shapeControlMesh_.cellCentres();

    Info<< "    Created cell centres" << endl;

    const PtrList<cellSizeAndAlignmentControl>& controlFunctions =
        sizeAndAlignment_.controlFunctions();

    DynamicList<Vb> verts(shapeControlMesh_.number_of_vertices());

    forAll(cellCentres, cellI)
    {
        const Foam::point& pt = cellCentres[cellI];

        if (!geometryToConformTo_.inside(pt))
        {
            continue;
        }

        scalarList bary;
        cellShapeControlMesh::Cell_handle ch;

        shapeControlMesh_.barycentricCoords(pt, bary, ch);

        if (shapeControlMesh_.is_infinite(ch))
        {
            continue;
        }

        scalar interpolatedSize = 0;
        forAll(bary, pI)
        {
            interpolatedSize += bary[pI]*ch->vertex(pI)->targetCellSize();
        }

        label lastPriority = labelMax;
        scalar lastCellSize = GREAT;
        forAll(controlFunctions, fI)
        {
            const cellSizeAndAlignmentControl& controlFunction =
                controlFunctions[fI];

            if (controlFunction.priority() > lastPriority)
            {
                continue;
            }

            if (isA<searchableSurfaceControl>(controlFunction))
            {
                const cellSizeFunction& sizeFunction =
                    dynamicCast<const searchableSurfaceControl>
                    (
                        controlFunction
                    ).sizeFunction();

                scalar cellSize = 0;
                if (sizeFunction.cellSize(pt, cellSize))
                {
                    if (controlFunction.priority() == lastPriority)
                    {
                        if (cellSize < lastCellSize)
                        {
                            lastCellSize = cellSize;
                        }
                    }
                    else
                    {
                        lastCellSize = cellSize;
                    }

                    lastPriority = controlFunction.priority();
                }
            }
        }

        if
        (
            lastCellSize < GREAT
         && mag(interpolatedSize - lastCellSize)/lastCellSize > 0.2
        )
        {
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
            verts.last().targetCellSize() = lastCellSize;
            verts.last().alignment() = tensor::I;
        }
    }

    shapeControlMesh_.insertPoints(verts);

    return verts.size();
}


void Foam::cellShapeControl::smoothMesh()
{
    label maxSmoothingIterations = 200;
    scalar minResidual = 0;

    labelListList pointPoints;
    autoPtr<mapDistribute> meshDistributor = buildMap
    (
        shapeControlMesh_,
        pointPoints
    );

    triadField alignments = buildAlignmentField(shapeControlMesh_);
    pointField points = buildPointField(shapeControlMesh_);
    // Setup the sizes and alignments on each point
    triadField fixedAlignments(shapeControlMesh_.vertexCount(), triad::unset);

    for
    (
        CellSizeDelaunay::Finite_vertices_iterator vit =
            shapeControlMesh_.finite_vertices_begin();
        vit != shapeControlMesh_.finite_vertices_end();
        ++vit
    )
    {
        if (vit->real())
        {
            const tensor& alignment = vit->alignment();

            fixedAlignments[vit->index()] = triad
            (
                alignment.x(),
                alignment.y(),
                alignment.z()
            );
        }
    }

    Info<< nl << "Smoothing alignments" << endl;

    for (label iter = 0; iter < maxSmoothingIterations; iter++)
    {
        Info<< "Iteration " << iter;

        meshDistributor().distribute(points);
        meshDistributor().distribute(alignments);

        scalar residual = 0;

        triadField triadAv(alignments.size(), triad::unset);

        forAll(pointPoints, pI)
        {
            const labelList& pPoints = pointPoints[pI];

            if (pPoints.empty())
            {
                continue;
            }

            triad& newTriad = triadAv[pI];

            forAll(pPoints, adjPointI)
            {
                const label adjPointIndex = pPoints[adjPointI];

                scalar dist = mag(points[pI] - points[adjPointIndex]);

                dist = max(dist, SMALL);

                triad tmpTriad = alignments[adjPointIndex];

                for (direction dir = 0; dir < 3; dir++)
                {
                    if (tmpTriad.set(dir))
                    {
                        tmpTriad[dir] *= (1.0/dist);
                    }
                }

                newTriad += tmpTriad;
            }

            newTriad.normalize();
            newTriad.orthogonalize();
            newTriad = newTriad.sortxyz();

            // Enforce the boundary conditions
            const triad& fixedAlignment = fixedAlignments[pI];

            label nFixed = 0;

            forAll(fixedAlignment, dirI)
            {
                if (fixedAlignment[dirI] != triad::unset[dirI])
                {
                    nFixed++;
                }
            }

            if (nFixed == 1)
            {
                forAll(fixedAlignment, dirI)
                {
                    if (fixedAlignment.set(dirI))
                    {
                        newTriad.align(fixedAlignment[dirI]);
                    }
                }
            }
            else if (nFixed == 2)
            {
                forAll(fixedAlignment, dirI)
                {
                    if (fixedAlignment.set(dirI))
                    {
                        newTriad[dirI] = fixedAlignment[dirI];
                    }
                    else
                    {
                        newTriad[dirI] = triad::unset[dirI];
                    }
                }

                newTriad.orthogonalize();
            }
            else if (nFixed == 3)
            {
                forAll(fixedAlignment, dirI)
                {
                    if (fixedAlignment.set(dirI))
                    {
                        newTriad[dirI] = fixedAlignment[dirI];
                    }
                }
            }

            const triad& oldTriad = alignments[pI];

            if (newTriad.set(vector::X) && oldTriad.set(vector::X))
            {
                scalar dotProd = (oldTriad.x() & newTriad.x());

                scalar diff = mag(dotProd) - 1.0;
                residual += mag(diff);
            }
            if (newTriad.set(vector::Y) && oldTriad.set(vector::Y))
            {
                scalar dotProd = (oldTriad.y() & newTriad.y());

                scalar diff = mag(dotProd) - 1.0;
                residual += mag(diff);
            }
            if (newTriad.set(vector::Z) && oldTriad.set(vector::Z))
            {
                scalar dotProd = (oldTriad.z() & newTriad.z());

                scalar diff = mag(dotProd) - 1.0;
                residual += mag(diff);
            }
        }

        forAll(alignments, pI)
        {
            alignments[pI] = triadAv[pI].sortxyz();
        }

        reduce(residual, sumOp<scalar>());

        Info<< ", Residual = " << residual << endl;

        if (residual <= minResidual)
        {
            break;
        }
    }

    meshDistributor().distribute(alignments);

    for
    (
        CellSizeDelaunay::Finite_vertices_iterator vit =
            shapeControlMesh_.finite_vertices_begin();
        vit != shapeControlMesh_.finite_vertices_end();
        ++vit
    )
    {
        if (vit->real())
        {
            vit->alignment() = tensor
            (
                alignments[vit->index()].x(),
                alignments[vit->index()].y(),
                alignments[vit->index()].z()
            );
        }
    }

    labelList referredPoints;
    autoPtr<mapDistribute> referredDistributor = buildReferredMap
    (
        shapeControlMesh_,
        referredPoints
    );

    alignments.setSize(shapeControlMesh_.vertexCount());
    referredDistributor().distribute(alignments);

    label referredI = 0;
    for
    (
        CellSizeDelaunay::Finite_vertices_iterator vit =
            shapeControlMesh_.finite_vertices_begin();
        vit != shapeControlMesh_.finite_vertices_end();
        ++vit
    )
    {
        if (vit->referred())
        {
            const triad& t = alignments[referredPoints[referredI++]];

            vit->alignment() = tensor(t.x(), t.y(), t.z());
        }
    }
}


// ************************************************************************* //
