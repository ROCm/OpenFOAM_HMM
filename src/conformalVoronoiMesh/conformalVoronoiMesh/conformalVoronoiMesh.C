/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "conformalVoronoiMesh.H"
#include "initialPointsMethod.H"
#include "uint.H"
#include "ulong.H"
#include "surfaceFeatures.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::conformalVoronoiMesh::dualCellSurfaceIntersection
(
    const Triangulation::Finite_vertices_iterator& vit
) const
{
    std::list<Facet>  facets;
    incident_facets(vit, std::back_inserter(facets));

    for
    (
        std::list<Facet>::iterator fit=facets.begin();
        fit != facets.end();
        ++fit
    )
    {
        if
        (
            is_infinite(fit->first)
            || is_infinite(fit->first->neighbor(fit->second))
        )
        {
            return true;
        }

        point dE0 = topoint(dual(fit->first));

        // If edge end is outside bounding box then edge cuts boundary
        if (!geometryToConformTo_.bounds().contains(dE0))
        {
            return true;
        }

        point dE1 = topoint(dual(fit->first->neighbor(fit->second)));

        // If other edge end is outside bounding box then edge cuts boundary
        if (!geometryToConformTo_.bounds().contains(dE1))
        {
            return true;
        }

        // Check for the edge passing through a surface
        if (geometryToConformTo_.findSurfaceAnyIntersection(dE0, dE1))
        {
            return true;
        }
    }

    return false;
}


void Foam::conformalVoronoiMesh::calcDualMesh
(
    pointField& points,
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    wordList& patchNames,
    labelList& patchSizes,
    labelList& patchStarts
)
{
    // ~~~~~~~~~~~ removing short edges by indexing dual vertices ~~~~~~~~~~~~~~

    for
    (
        Triangulation::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        cit->cellIndex() = -1;
    }

    points.setSize(number_of_cells());

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Looking up details from a dictionary, in future the will be available
    // from a controls class.

    const dictionary& cvMeshDict( cvMeshControls_.cvMeshDict());

    scalar defaultCellSize
    (
        readScalar
        (
            cvMeshDict.subDict("motionControl").lookup("defaultCellSize")
        )
    );

    scalar minimumEdgeLengthCoeff
    (
        readScalar
        (
            cvMeshDict.subDict("polyMeshFiltering").lookup
            (
                "minimumEdgeLengthCoeff"
            )
        )
    );

    scalar minEdgeLenSqr = sqr(defaultCellSize*minimumEdgeLengthCoeff);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    label dualVerti = 0;

    // Scanning by number of short (dual) edges (nSE) attached to the
    // circumcentre of each Delaunay tet.  A Delaunay tet may only have four
    // dual edges emanating from its circumcentre, assigning positions and
    // indices to those with 4 short edges attached first, then >= 3, then >= 2
    // etc.
    for (label nSE = 4; nSE >= 0; nSE--)
    {
        Info<< nl << "Scanning for dual vertices with >= "
            << nSE
            << " short edges attached." << endl;

        for
        (
            Triangulation::Finite_cells_iterator cit = finite_cells_begin();
            cit != finite_cells_end();
            ++cit
        )
        {
            // If the Delaunay tet has an index already then it has either
            // evaluated itself and taken action or has had its index dictated
            // by a neighbouring tet with more short edges attached.

            if (cit->cellIndex() == -1)
            {
                point dualVertex = topoint(dual(cit));

                label shortEdges = 0;

                List<bool> edgeIsShort(4, false);

                List<bool> neighbourAlreadyIndexed(4, false);

                // Loop over the four facets of the Delaunay tet
                for (label f = 0; f < 4; f++)
                {
                    // Check that at least one of the vertices of the facet is
                    // an internal or boundary point
                    if
                    (
                        cit->vertex(vertex_triple_index(f, 0))->
                        internalOrBoundaryPoint()
                        || cit->vertex(vertex_triple_index(f, 1))->
                        internalOrBoundaryPoint()
                        || cit->vertex(vertex_triple_index(f, 2))->
                        internalOrBoundaryPoint()
                    )
                    {
                        point neighDualVertex;

                        label cNI = cit->neighbor(f)->cellIndex();

                        if (cNI == -1)
                        {
                            neighDualVertex = topoint(dual(cit->neighbor(f)));
                        }
                        else
                        {
                            neighDualVertex = points[cNI];
                        }

                        if
                        (
                            magSqr(dualVertex - neighDualVertex) < minEdgeLenSqr
                        )
                        {
                            edgeIsShort[f] = true;

                            if (cNI > -1)
                            {
                                neighbourAlreadyIndexed[f] = true;
                            }

                            shortEdges++;
                        }
                    }
                }

                if (nSE == 0 && shortEdges == 0)
                {
                    // Final iteration and no short edges are found, index
                    // remaining dual vertices.

                    if
                    (
                        cit->vertex(0)->internalOrBoundaryPoint()
                     || cit->vertex(1)->internalOrBoundaryPoint()
                     || cit->vertex(2)->internalOrBoundaryPoint()
                     || cit->vertex(3)->internalOrBoundaryPoint()
                    )
                    {
                        cit->cellIndex() = dualVerti;
                        points[dualVerti] = dualVertex;
                        dualVerti++;
                    }
                }
                else if
                (
                    shortEdges >= nSE
                )
                {
                    // Info<< neighbourAlreadyIndexed << ' '
                    //     << edgeIsShort << endl;

                    label numUnindexedNeighbours = 1;

                    for (label f = 0; f < 4; f++)
                    {
                        if (edgeIsShort[f] && !neighbourAlreadyIndexed[f])
                        {
                            dualVertex += topoint(dual(cit->neighbor(f)));

                            numUnindexedNeighbours++;
                        }
                    }

                    dualVertex /= numUnindexedNeighbours;

                    label nearestExistingIndex = -1;

                    point nearestIndexedNeighbourPos = vector::zero;

                    scalar minDistSqrToNearestIndexedNeighbour = VGREAT;

                    for (label f = 0; f < 4; f++)
                    {
                        if (edgeIsShort[f] && neighbourAlreadyIndexed[f])
                        {
                            label cNI = cit->neighbor(f)->cellIndex();

                            point indexedNeighbourPos = points[cNI];

                            if
                            (
                                magSqr(indexedNeighbourPos - dualVertex)
                              < minDistSqrToNearestIndexedNeighbour
                            )
                            {
                                nearestExistingIndex = cNI;

                                nearestIndexedNeighbourPos =
                                indexedNeighbourPos;

                                minDistSqrToNearestIndexedNeighbour =
                                magSqr(indexedNeighbourPos - dualVertex);
                            }
                        }
                    }

                    if
                    (
                        nearestExistingIndex > -1
                     && minDistSqrToNearestIndexedNeighbour < minEdgeLenSqr
                    )
                    {
                        points[nearestExistingIndex] =
                        0.5*(dualVertex + nearestIndexedNeighbourPos);

                        for (label f = 0; f < 4; f++)
                        {
                            if (edgeIsShort[f] && !neighbourAlreadyIndexed[f])
                            {
                                cit->neighbor(f)->cellIndex() =
                                nearestExistingIndex;
                            }
                        }

                        cit->cellIndex() = nearestExistingIndex;
                    }
                    else
                    {
                        for (label f = 0; f < 4; f++)
                        {
                            if (edgeIsShort[f] && !neighbourAlreadyIndexed[f])
                            {
                                cit->neighbor(f)->cellIndex() = dualVerti;
                            }
                        }

                        cit->cellIndex() = dualVerti;

                        points[dualVerti] = dualVertex;

                        dualVerti++;
                    }
                }
            }
        }
    }

    points.setSize(dualVerti);

    // ~~~~~~~~~~~~~~~~~~~~~~~~~ dual cell indexing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // assigns an index to the Delaunay vertices which will be the dual cell
    // index used for owner neighbour assignment.

    // The indices of the points are reset which destroys the point-pair
    // matching, so the type of each vertex are reset to avoid any ambiguity.

    label dualCelli = 0;

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalOrBoundaryPoint())
        {
            vit->type() = Vb::INTERNAL_POINT;
            vit->index() = dualCelli;
            dualCelli++;
        }
        else
        {
            vit->type() = Vb::FAR_POINT;
            vit->index() = -1;
        }
    }

    // ~~~~~~~~~~~~ dual face and owner neighbour construction ~~~~~~~~~~~~~~~~~

    //label nPatches = qSurf_.patches().size() + 1;

    //label defaultPatchIndex = qSurf_.patches().size();

    label nPatches = 1;

    label defaultPatchIndex = 0;

    patchNames.setSize(nPatches);

    //const geometricSurfacePatchList& surfacePatches = qSurf_.patches();

    // forAll(surfacePatches, sP)
    // {
    //     patchNames[sP] = surfacePatches[sP].name();
    // }

    patchNames[defaultPatchIndex] = "cvMesh_defaultPatch";

    patchSizes.setSize(nPatches);

    patchStarts.setSize(nPatches);

    List<DynamicList<face> > patchFaces(nPatches, DynamicList<face>(0));

    List<DynamicList<label> > patchOwners(nPatches, DynamicList<label>(0));

    faces.setSize(number_of_edges());

    owner.setSize(number_of_edges());

    neighbour.setSize(number_of_edges());

    label dualFacei = 0;

    for
    (
        Triangulation::Finite_edges_iterator eit = finite_edges_begin();
        eit != finite_edges_end();
        ++eit
    )
    {
        Cell_handle c = eit->first;
        Vertex_handle vA = c->vertex(eit->second);
        Vertex_handle vB = c->vertex(eit->third);

        if
        (
            vA->internalOrBoundaryPoint()
         || vB->internalOrBoundaryPoint()
        )
        {
            Cell_circulator ccStart = incident_cells(*eit);
            Cell_circulator cc1 = ccStart;
            Cell_circulator cc2 = cc1;

            // Advance the second circulator so that it always stays on the next
            // cell around the edge;
            cc2++;

            DynamicList<label> verticesOnFace;

            do
            {
                label cc1I = cc1->cellIndex();

                label cc2I = cc2->cellIndex();


                if (cc1I < 0 || cc2I < 0)
                {
                    FatalErrorIn("Foam::conformalVoronoiMesh::calcDualMesh")
                        << "Dual face uses circumcenter defined by a "
                        << "Delaunay tetrahedron with no internal "
                        << "or boundary points.  Defining Delaunay edge ends: "
                        << topoint(vA->point()) << " "
                        << topoint(vB->point()) << nl
                        << exit(FatalError);
                }

                if (cc1I != cc2I)
                {
                    verticesOnFace.append(cc1I);
                }

                cc1++;

                cc2++;
            } while (cc1 != ccStart);

            verticesOnFace.shrink();

            if (verticesOnFace.size() >= 3)
            {
                face newDualFace(verticesOnFace);

                label dcA = vA->index();

                if (!vA->internalOrBoundaryPoint())
                {
                    dcA = -1;
                }

                label dcB = vB->index();

                if (!vB->internalOrBoundaryPoint())
                {
                    dcB = -1;
                }

                label dcOwn = -1;
                label dcNei = -1;

                if (dcA == -1 && dcB == -1)
                {
                    FatalErrorIn("calcDualMesh")
                        << "Attempting to create a face joining "
                        << "two external dual cells "
                        << exit(FatalError);
                }
                else if (dcA == -1 || dcB == -1)
                {
                    // boundary face, find which is the owner

                    if (dcA == -1)
                    {
                        dcOwn = dcB;

                        // reverse face order to correctly orientate normal
                        reverse(newDualFace);
                    }
                    else
                    {
                        dcOwn = dcA;
                    }

                    // Find which patch this face is on by finding the
                    // intersection with the surface of the Delaunay edge
                    // generating the face and identify the region of the
                    // intersection.

                    point ptA = topoint(vA->point());

                    point ptB = topoint(vB->point());

                    //pointIndexHit pHit = qSurf_.tree().findLineAny(ptA, ptB);

                    //label patchIndex = qSurf_[pHit.index()].region();

                    label patchIndex = defaultPatchIndex;

                    if (patchIndex == -1)
                    {
                        patchIndex = defaultPatchIndex;

                        WarningIn("Foam::conformalVoronoiMesh::calcDualMesh")
                            << "Dual face found that is not on a surface "
                            << "patch. Adding to "
                            << patchNames[defaultPatchIndex]
                            << endl;
                    }

                    patchFaces[patchIndex].append(newDualFace);
                    patchOwners[patchIndex].append(dcOwn);
                }
                else
                {
                    // internal face, find the lower cell to be the owner

                    if (dcB > dcA)
                    {
                        dcOwn = dcA;
                        dcNei = dcB;
                    }
                    else
                    {
                        dcOwn = dcB;
                        dcNei = dcA;

                        // reverse face order to correctly orientate normal
                        reverse(newDualFace);
                    }

                    faces[dualFacei] = newDualFace;

                    owner[dualFacei] = dcOwn;

                    neighbour[dualFacei] = dcNei;

                    dualFacei++;
                }
            }
            // else
            // {
            //     Info<< verticesOnFace.size()
            //         << " size face not created." << endl;
            // }
        }
    }

    label nInternalFaces = dualFacei;

    faces.setSize(nInternalFaces);

    owner.setSize(nInternalFaces);

    neighbour.setSize(nInternalFaces);

    // ~~~~~~~~ sort owner, reordinging neighbour and faces to match ~~~~~~~~~~~
    // two stage sort for upper triangular order:  sort by owner first, then for
    // each block of owners sort by neighbour

    labelList sortingIndices;

    // Stage 1

    {
        SortableList<label> sortedOwner(owner);

        sortingIndices = sortedOwner.indices();
    }

    {
        labelList copyOwner(owner.size());

        forAll(sortingIndices, sI)
        {
            copyOwner[sI] = owner[sortingIndices[sI]];
        }

        owner = copyOwner;
    }

    {
        labelList copyNeighbour(neighbour.size());

        forAll(sortingIndices, sI)
        {
            copyNeighbour[sI] = neighbour[sortingIndices[sI]];
        }

        neighbour = copyNeighbour;
    }

    {
        faceList copyFaces(faces.size());

        forAll(sortingIndices, sI)
        {
            copyFaces[sI] = faces[sortingIndices[sI]];
        }

        faces = copyFaces;
    }

    // Stage 2

    sortingIndices = -1;

    DynamicList<label> ownerCellJumps;

    // Force first owner entry to be a jump
    ownerCellJumps.append(0);

    for (label o = 1; o < owner.size(); o++)
    {
        if (owner[o] > owner[o-1])
        {
            ownerCellJumps.append(o);
        }
    }

    ownerCellJumps.shrink();

    forAll(ownerCellJumps, oCJ)
    {
        label start = ownerCellJumps[oCJ];

        label length;

        if (oCJ == ownerCellJumps.size() - 1)
        {
            length = owner.size() - start;
        }
        else
        {
            length = ownerCellJumps[oCJ + 1] - start;
        }

        SubList<label> neighbourBlock(neighbour, length, start);

        SortableList<label> sortedNeighbourBlock(neighbourBlock);

        forAll(sortedNeighbourBlock, sNB)
        {
            sortingIndices[start + sNB] =
            sortedNeighbourBlock.indices()[sNB] + start;
        }
    }

    // Perform sort

    {
        labelList copyOwner(owner.size());

        forAll(sortingIndices, sI)
        {
            copyOwner[sI] = owner[sortingIndices[sI]];
        }

        owner = copyOwner;
    }

    {
        labelList copyNeighbour(neighbour.size());

        forAll(sortingIndices, sI)
        {
            copyNeighbour[sI] = neighbour[sortingIndices[sI]];
        }

        neighbour = copyNeighbour;
    }

    {
        faceList copyFaces(faces.size());

        forAll(sortingIndices, sI)
        {
            copyFaces[sI] = faces[sortingIndices[sI]];
        }

        faces = copyFaces;
    }

    // ~~~~~~~~ add patch information ~~~~~~~~~~~

    label nBoundaryFaces = 0;

    forAll(patchFaces, p)
    {
        patchFaces[p].shrink();

        patchOwners[p].shrink();

        patchSizes[p] = patchFaces[p].size();

        patchStarts[p] = nInternalFaces + nBoundaryFaces;

        nBoundaryFaces += patchSizes[p];
    }

    faces.setSize(nInternalFaces + nBoundaryFaces);

    owner.setSize(nInternalFaces + nBoundaryFaces);

    forAll(patchFaces, p)
    {
        forAll(patchFaces[p], f)
        {
            faces[dualFacei] = patchFaces[p][f];

            owner[dualFacei] = patchOwners[p][f];

            dualFacei++;
        }
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conformalVoronoiMesh::conformalVoronoiMesh
(
    const Time& runTime,
    const IOdictionary& cvMeshDict
)
:
    HTriangulation(),
    runTime_(runTime),
    allGeometry_
    (
        IOobject
        (
            "cvSearchableSurfacesDirectory",
            runTime_.constant(),
            "triSurface",
            runTime_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        cvMeshDict.subDict("geometry")
    ),
    geometryToConformTo_
    (
        *this,
        allGeometry_,
        cvMeshDict.subDict("surfaceConformation")
    ),
    cvMeshControls_(*this, cvMeshDict),
    startOfInternalPoints_(0),
    startOfSurfacePointPairs_(0),
    initialPointsMethod_
    (
        initialPointsMethod::New
        (
            cvMeshDict.subDict("initialPoints"),
            *this
        )
    )
{
    timeCheck();

    conformToFeaturePoints();
    timeCheck();

    insertInitialPoints();
    timeCheck();

    conformToSurface();
    timeCheck();

    writePoints("allPoints.obj", false);
    timeCheck();

    writeMesh();
    timeCheck();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conformalVoronoiMesh::~conformalVoronoiMesh()
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::conformalVoronoiMesh::timeCheck() const
{
    Info<< nl << "--- [ " << runTime_.elapsedCpuTime() << "s, delta "
        << runTime_.cpuTimeIncrement()<< "s ] --- " << endl;
}


void Foam::conformalVoronoiMesh::insertSurfacePointPairs
(
    const List<scalar>& surfacePpDist,
    const List<point>& surfacePoints,
    const List<vector>& surfaceNormals,
    const fileName fName
)
{
    if
    (
        surfacePpDist.size() != surfacePoints.size()
     || surfacePpDist.size() != surfaceNormals.size()
    )
    {
        FatalErrorIn("Foam::conformalVoronoiMesh::insertPointPairs")
            << "surfacePpDist, surfacePoints and surfaceNormals are not "
            << "the same size. Sizes"
            << surfacePpDist.size() << ' '
            << surfacePoints.size() << ' '
            << surfaceNormals.size()
            << exit(FatalError);
    }

    forAll(surfacePoints, p)
    {
        insertPointPair
        (
            surfacePpDist[p],
            surfacePoints[p],
            surfaceNormals[p]
        );
    }

    if (fName != fileName::null)
    {
        writePoints(fName, surfacePoints);
    }
}


void Foam::conformalVoronoiMesh::conformToFeaturePoints()
{
    Info<< nl << "Conforming to feature points" << endl;



    Info<< "   Conforming to " << "XXX" << " feature locations" << nl
        << "   Inserting " << "YYY" << " points" << endl;
}


void Foam::conformalVoronoiMesh::reinsertFeaturePoints()
{

}


void Foam::conformalVoronoiMesh::insertInitialPoints()
{
    startOfInternalPoints_ = number_of_vertices();

    label nVert = startOfInternalPoints_;

    Info<< nl << "Inserting initial points" << endl;

    std::vector<Point> initialPoints = initialPointsMethod_->initialPoints();

    Info<< "    " << initialPoints.size() << " points to insert..." << endl;

    // using the range insert (faster than inserting points one by one)
    insert(initialPoints.begin(), initialPoints.end());

    Info<< "    " << number_of_vertices() - startOfInternalPoints_
        << " points inserted" << endl;

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->uninitialised())
        {
            vit->index() = nVert++;
        }
    }

    writePoints("initialPoints.obj", true);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::conformalVoronoiMesh::conformToSurface()
{
    Info<< nl << "Conforming to surfaces" << endl;

    startOfSurfacePointPairs_ = number_of_vertices();

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Looking up details from a dictionary, in future the will be available
    // from a controls class.

    const dictionary& cvMeshDict( cvMeshControls_.cvMeshDict());

    scalar defaultCellSize
    (
        readScalar
        (
            cvMeshDict.subDict("motionControl").lookup("defaultCellSize")
        )
    );

    scalar surfDepthCoeff
    (
        readScalar
        (
            cvMeshDict.subDict("surfaceConformation").lookup
            (
                "surfacePointSearchDepthCoeff"
            )
        )
    );

    scalar ppDistCoeff
    (
        readScalar
        (
            cvMeshDict.subDict("surfaceConformation").lookup
            (
                "pointPairDistanceCoeff"
            )
        )
    );

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Surface protrusion conformation

    label nIterations = 1;

    for(label iterationNo = 0; iterationNo < nIterations; iterationNo++)
    {
        DynamicList<scalar> surfacePpDist;
        DynamicList<point> surfacePoints;
        DynamicList<vector> surfaceNormals;

        for
        (
            Triangulation::Finite_vertices_iterator vit =
            finite_vertices_begin();
            vit != finite_vertices_end();
            vit++
        )
        {
            if (vit->internalPoint())
            {
                point vert(topoint(vit->point()));

                // TODO Need to have a function to recover the local cell size,
                // use the defaultCellSize for the moment

                scalar searchDistanceSqr = sqr(defaultCellSize*surfDepthCoeff);
                pointIndexHit pHit;
                vector normal;

                geometryToConformTo_.findSurfaceNearestAndNormal
                (
                    vert,
                    searchDistanceSqr,
                    pHit,
                    normal
                );

                if (pHit.hit())
                {
                    vit->setNearBoundary();

                    if (dualCellSurfaceIntersection(vit))
                    {
                        // If the point is within a given distance of a feature
                        // edge, shift it to being an edge control point
                        // instead, this will prevent "pits" forming.

                        surfacePpDist.append(defaultCellSize*ppDistCoeff);
                        surfacePoints.append(pHit.hitPoint());
                        surfaceNormals.append(normal);
                    }
                }
            }
        }

        Info<< nl <<"    iterationNo " << iterationNo << nl
            << "    number_of_vertices " << number_of_vertices() << nl
            << "    surfacePoints.size() " << surfacePoints.size() << endl;

        insertSurfacePointPairs
        (
            surfacePpDist,
            surfacePoints,
            surfaceNormals,
            fileName
            (
                "surfaceConformationLocations_" + name(iterationNo) + ".obj"
            )
        );

        // Feature edge conformation.

        geometryToConformTo_.findEdgeNearest
        (
            pointField(surfacePoints),
            scalarField
            (
                surfacePoints.size(),
                sqr(defaultCellSize*surfDepthCoeff)
            )
        );

        // After the surface conformation points are added, any points that are
        // still protruding the surface may be protruding from edges, so
        // identify the points and test if they are close to a feature edge

        DynamicList<point> edgeGenerationPoints;

        for
        (
            Triangulation::Finite_vertices_iterator vit =
            finite_vertices_begin();
            vit != finite_vertices_end();
            vit++
        )
        {
            if (vit->nearBoundary())
            {
                if (dualCellSurfaceIntersection(vit))
                {
                    // Test to see if near to an edge, conform to the nearest
                    // point on that edge if so.

                    edgeGenerationPoints.append(topoint(vit->point()));
                }
            }
        }

        writePoints
        (
            "edgeGenerationLocations_" + name(iterationNo) + ".obj",
            edgeGenerationPoints
        );

    }
}


// ************************************************************************* //
