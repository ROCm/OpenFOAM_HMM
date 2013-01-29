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
#include "IOstreams.H"
#include "OFstream.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "ListOps.H"
#include "polyMeshFilter.H"
#include "polyTopoChange.H"
#include "PrintTable.H"
#include "pointMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::conformalVoronoiMesh::timeCheck
(
    const string& description
) const
{
    if (cvMeshControls().timeChecks())
    {
        Info<< nl << "--- [ cpuTime "
            << runTime_.elapsedCpuTime() << " s, "
            << "delta " << runTime_.cpuTimeIncrement()<< " s";

        if (description != word::null)
        {
            Info<< ", " << description << " ";
        }
        else
        {
            Info<< " ";
        }

        Info<< "] --- " << endl;

        memInfo m;

        if (m.valid())
        {
            PrintTable<word, label> memoryTable("Memory Usage (kB)");

            memoryTable.add("mSize", m.size());
            memoryTable.add("mPeak", m.peak());
            memoryTable.add("mRss", m.rss());

            Info<< incrIndent;
            memoryTable.print(Info);
            Info<< decrIndent;
        }
    }
}


void Foam::conformalVoronoiMesh::printVertexInfo() const
{
    label nInternal = 0;
    label nInternalRef = 0;
    label nUnassigned = 0;
    label nUnassignedRef = 0;
    label nInternalNearBoundary = 0;
    label nInternalNearBoundaryRef = 0;
    label nInternalSurface = 0;
    label nInternalSurfaceRef = 0;
    label nInternalFeatureEdge = 0;
    label nInternalFeatureEdgeRef = 0;
    label nInternalFeaturePoint = 0;
    label nInternalFeaturePointRef = 0;
    label nExternalSurface = 0;
    label nExternalSurfaceRef = 0;
    label nExternalFeatureEdge = 0;
    label nExternalFeatureEdgeRef = 0;
    label nExternalFeaturePoint = 0;
    label nExternalFeaturePointRef = 0;
    label nFar = 0;
    label nReferred = 0;

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->type() == Vb::vtInternal)
        {
            if (vit->referred())
            {
                nReferred++;
                nInternalRef++;
            }

            nInternal++;
        }
        else if (vit->type() == Vb::vtUnassigned)
        {
            if (vit->referred())
            {
                nReferred++;
                nUnassignedRef++;
            }

            nUnassigned++;
        }
        else if (vit->type() == Vb::vtInternalNearBoundary)
        {
            if (vit->referred())
            {
                nReferred++;
                nInternalNearBoundaryRef++;
            }

            nInternalNearBoundary++;
        }
        else if (vit->type() == Vb::vtInternalSurface)
        {
            if (vit->referred())
            {
                nReferred++;
                nInternalSurfaceRef++;
            }

            nInternalSurface++;
        }
        else if (vit->type() == Vb::vtInternalFeatureEdge)
        {
            if (vit->referred())
            {
                nReferred++;
                nInternalFeatureEdgeRef++;
            }

            nInternalFeatureEdge++;
        }
        else if (vit->type() == Vb::vtInternalFeaturePoint)
        {
            if (vit->referred())
            {
                nReferred++;
                nInternalFeaturePointRef++;
            }

            nInternalFeaturePoint++;
        }
        else if (vit->type() == Vb::vtExternalSurface)
        {
            if (vit->referred())
            {
                nReferred++;
                nExternalSurfaceRef++;
            }

            nExternalSurface++;
        }
        else if (vit->type() == Vb::vtExternalFeatureEdge)
        {
            if (vit->referred())
            {
                nReferred++;
                nExternalFeatureEdgeRef++;
            }

            nExternalFeatureEdge++;
        }
        else if (vit->type() == Vb::vtExternalFeaturePoint)
        {
            if (vit->referred())
            {
                nReferred++;
                nExternalFeaturePointRef++;
            }

            nExternalFeaturePoint++;
        }
        else if (vit->type() == Vb::vtFar)
        {
            nFar++;
        }
    }

    label nTotalVertices
        = nUnassigned
        + nInternal
        + nInternalNearBoundary
        + nInternalSurface
        + nInternalFeatureEdge
        + nInternalFeaturePoint
        + nExternalSurface
        + nExternalFeatureEdge
        + nExternalFeaturePoint
        + nFar;

    if (nTotalVertices != label(number_of_vertices()))
    {
        WarningIn("Foam::conformalVoronoiMesh::printVertexInfo()")
            << nTotalVertices << " does not equal " << number_of_vertices()
            << endl;
    }

    PrintTable<word, label> vertexTable("Vertex Type Information");

    vertexTable.add("Total", nTotalVertices);
    vertexTable.add("Unassigned", nUnassigned);
    vertexTable.add("nInternal", nInternal);
    vertexTable.add("nInternalNearBoundary", nInternalNearBoundary);
    vertexTable.add("nInternalSurface", nInternalSurface);
    vertexTable.add("nInternalFeatureEdge", nInternalFeatureEdge);
    vertexTable.add("nInternalFeaturePoint", nInternalFeaturePoint);
    vertexTable.add("nExternalSurface", nExternalSurface);
    vertexTable.add("nExternalFeatureEdge", nExternalFeatureEdge);
    vertexTable.add("nExternalFeaturePoint", nExternalFeaturePoint);
    vertexTable.add("nFar", nFar);
    vertexTable.add("nReferred", nReferred);

    Info<< endl;
    vertexTable.print(Info);
}


void Foam::conformalVoronoiMesh::drawDelaunayCell
(
    Ostream& os,
    const Cell_handle& c,
    label offset
) const
{
    // Supply offset as tet number
    offset *= 4;

    os  << "# cell index: " << label(c->cellIndex()) << endl;

    os  << "# circumradius "
        << mag(c->dual() - topoint(c->vertex(0)->point()))
        << endl;

    for (int i = 0; i < 4; i++)
    {
        os  << "# index / type / procIndex: "
            << label(c->vertex(i)->index()) << " "
            << label(c->vertex(i)->type()) << " "
            << label(c->vertex(i)->procIndex()) << endl;

        meshTools::writeOBJ(os, topoint(c->vertex(i)->point()));
    }

    os  << "f " << 1 + offset << " " << 3 + offset << " " << 2 + offset << nl
        << "f " << 2 + offset << " " << 3 + offset << " " << 4 + offset << nl
        << "f " << 1 + offset << " " << 4 + offset << " " << 3 + offset << nl
        << "f " << 1 + offset << " " << 2 + offset << " " << 4 + offset << endl;

//    os  << "# cicumcentre " << endl;

//    meshTools::writeOBJ(os, c->dual());

//    os  << "l " << 1 + offset << " " << 5 + offset << endl;
}


void Foam::conformalVoronoiMesh::writePoints
(
    const fileName& fName,
    const Foam::indexedVertexEnum::vertexType startPointType,
    const Foam::indexedVertexEnum::vertexType endPointType
) const
{
    OFstream str(runTime_.path()/fName);

    Pout<< nl << "Writing points of types:" << nl;

    forAllConstIter
    (
        HashTable<int>,
        Foam::indexedVertexEnum::vertexTypeNames_,
        iter
    )
    {
        if (iter() >= startPointType && iter() <= endPointType)
        {
            Pout<< "    " << iter.key() << nl;
        }
    }

    Pout<< "to " << str.name() << endl;

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->type() >= startPointType && vit->type() <= endPointType)
        {
            meshTools::writeOBJ(str, topoint(vit->point()));
        }
    }
}


void Foam::conformalVoronoiMesh::writePoints
(
    const fileName& fName,
    const Foam::indexedVertexEnum::vertexType pointType
) const
{
    writePoints(fName, pointType, pointType);
}


void Foam::conformalVoronoiMesh::writeBoundaryPoints
(
    const fileName& fName
) const
{
    OFstream str(runTime_.path()/fName);

    Pout<< nl << "Writing boundary points to " << str.name() << endl;

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (!vit->internalPoint())
        {
            meshTools::writeOBJ(str, topoint(vit->point()));
        }
    }
}


void Foam::conformalVoronoiMesh::writePoints
(
    const fileName& fName,
    const List<Foam::point>& points
) const
{
    if (points.size())
    {
        OFstream str(runTime_.path()/fName);

        Pout<< nl << "Writing " << points.size() << " points from pointList to "
            << str.name() << endl;

        forAll(points, p)
        {
            meshTools::writeOBJ(str, points[p]);
        }
    }
}


void Foam::conformalVoronoiMesh::writePoints
(
    const fileName& fName,
    const List<Vb>& points
) const
{
    if (points.size())
    {
        OFstream str(runTime_.path()/fName);

        Pout<< nl << "Writing " << points.size() << " points from pointList to "
            << str.name() << endl;

        forAll(points, p)
        {
            meshTools::writeOBJ(str, topoint(points[p].point()));
        }
    }
}


void Foam::conformalVoronoiMesh::writeProcessorInterface
(
    const fileName& fName,
    const faceList& faces
) const
{
    OFstream str(runTime_.path()/fName);

    pointField points(number_of_finite_cells(), point::max);

    for
    (
        Delaunay::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        if (!cit->hasFarPoint() && !is_infinite(cit))
        {
            points[cit->cellIndex()] = cit->dual();
        }
    }

    meshTools::writeOBJ(str, faces, points);
}


void Foam::conformalVoronoiMesh::writeInternalDelaunayVertices
(
    const fileName& instance
) const
{
    pointField internalDelaunayVertices(number_of_vertices());

    label vertI = 0;

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalPoint())
        {
            internalDelaunayVertices[vertI++] = topoint(vit->point());
        }
    }

    internalDelaunayVertices.setSize(vertI);

    pointIOField internalDVs
    (
        IOobject
        (
            "internalDelaunayVertices",
            instance,
            runTime_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        internalDelaunayVertices
    );

    Info<< nl
        << "Writing " << internalDVs.name()
        << " to " << internalDVs.instance()
        << endl;

    internalDVs.write();
}


void Foam::conformalVoronoiMesh::writeMesh(const fileName& instance)
{
    writeInternalDelaunayVertices(instance);

    // Per cell the Delaunay vertex
    labelList cellToDelaunayVertex;
    // Per patch, per face the Delaunay vertex
    labelListList patchToDelaunayVertex;
    // Per patch the start of the dual faces
    labelList dualPatchStarts;

    {
        pointField points;
        labelList boundaryPts(number_of_finite_cells(), -1);
        faceList faces;
        labelList owner;
        labelList neighbour;
        wordList patchTypes;
        wordList patchNames;
        labelList patchSizes;
        labelList procNeighbours;
        pointField cellCentres;

        PackedBoolList boundaryFacesToRemove;

        calcDualMesh
        (
            points,
            boundaryPts,
            faces,
            owner,
            neighbour,
            patchTypes,
            patchNames,
            patchSizes,
            dualPatchStarts,
            procNeighbours,
            cellCentres,
            cellToDelaunayVertex,
            patchToDelaunayVertex,
            boundaryFacesToRemove
        );

        Info<< nl << "Writing polyMesh to " << instance << endl;

        writeMesh
        (
            Foam::polyMesh::defaultRegion,
            instance,
            points,
            boundaryPts,
            faces,
            owner,
            neighbour,
            patchTypes,
            patchNames,
            patchSizes,
            dualPatchStarts,
            procNeighbours,
            cellCentres,
            boundaryFacesToRemove
        );
    }

    if (cvMeshControls().writeTetDualMesh())
    {
        // Determine map from Delaunay vertex to Dual mesh
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // From all Delaunay vertices to cell (positive index)
        // or patch face (negative index)
        labelList vertexToDualAddressing(number_of_vertices(), 0);

        forAll(cellToDelaunayVertex, cellI)
        {
            label vertI = cellToDelaunayVertex[cellI];

            if (vertexToDualAddressing[vertI] != 0)
            {
                FatalErrorIn("conformalVoronoiMesh::writeMesh(..)")
                    << "Delaunay vertex " << vertI
                    << " from cell " << cellI
                    << " is already mapped to "
                    << vertexToDualAddressing[vertI]
                    << exit(FatalError);
            }
            vertexToDualAddressing[vertI] = cellI+1;
        }

        forAll(patchToDelaunayVertex, patchI)
        {
            const labelList& patchVertices = patchToDelaunayVertex[patchI];

            forAll(patchVertices, i)
            {
                label vertI = patchVertices[i];

                if (vertexToDualAddressing[vertI] > 0)
                {
                    FatalErrorIn("conformalVoronoiMesh::writeMesh(..)")
                        << "Delaunay vertex " << vertI
                        << " from patch " << patchI
                        << " local index " << i
                        << " is already mapped to cell "
                        << vertexToDualAddressing[vertI]-1
                        << exit(FatalError);
                }

                // Vertex might be used by multiple faces. Which one to
                // use? For now last one wins.
                label dualFaceI = dualPatchStarts[patchI]+i;
                vertexToDualAddressing[vertI] = -dualFaceI-1;
            }
        }


        // Calculate tet mesh addressing
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        pointField points;
        labelList boundaryPts(number_of_finite_cells(), -1);
        // From tet point back to Delaunay vertex index
        labelList pointToDelaunayVertex;
        faceList faces;
        labelList owner;
        labelList neighbour;
        wordList patchTypes;
        wordList patchNames;
        labelList patchSizes;
        labelList patchStarts;
        pointField cellCentres;

        calcTetMesh
        (
            points,
            pointToDelaunayVertex,
            faces,
            owner,
            neighbour,
            patchTypes,
            patchNames,
            patchSizes,
            patchStarts
        );



        // Calculate map from tet points to dual mesh cells/patch faces
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        labelIOList pointDualAddressing
        (
            IOobject
            (
                "pointDualAddressing",
                instance,
                "tetDualMesh"/polyMesh::meshSubDir,
                runTime_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            UIndirectList<label>
            (
                vertexToDualAddressing,
                pointToDelaunayVertex
            )()
        );

        label pointI = findIndex(pointDualAddressing, -1);
        if (pointI != -1)
        {
            WarningIn
            (
                "conformalVoronoiMesh::writeMesh\n"
                "(\n"
                "    const fileName& instance,\n"
                "    bool filterFaces\n"
                ")\n"
            )   << "Delaunay vertex " << pointI
                << " does not have a corresponding dual cell." << endl;
        }

        Info<< "Writing map from tetDualMesh points to Voronoi mesh to "
            << pointDualAddressing.objectPath() << endl;
        pointDualAddressing.write();



        // Write tet points corresponding to the Voronoi cell/face centre
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        {
            // Read Voronoi mesh
            fvMesh mesh
            (
                IOobject
                (
                    Foam::polyMesh::defaultRegion,
                    instance,
                    runTime_,
                    IOobject::MUST_READ
                )
            );
            pointIOField dualPoints
            (
                IOobject
                (
                    "dualPoints",
                    instance,
                    "tetDualMesh"/polyMesh::meshSubDir,
                    runTime_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false
                ),
                points
            );

            forAll(pointDualAddressing, pointI)
            {
                label index = pointDualAddressing[pointI];

                if (index > 0)
                {
                    label cellI = index-1;
                    dualPoints[pointI] = mesh.cellCentres()[cellI];
                }
                else if (index < 0)
                {
                    label faceI = -index-1;
                    if (faceI >= mesh.nInternalFaces())
                    {
                        dualPoints[pointI] = mesh.faceCentres()[faceI];
                    }
                }
            }

            Info<< "Writing new tetDualMesh points mapped onto Voronoi mesh to "
                << dualPoints.objectPath() << endl
                << "Replace the polyMesh/points with these." << endl;
            dualPoints.write();
        }


        labelList procNeighbours(patchNames.size(), -1);

//        Info<< nl << "Writing tetDualMesh to " << instance << endl;

//        writeMesh
//        (
//            "tetDualMesh",
//            instance,
//            points,
//            boundaryPts,
//            faces,
//            owner,
//            neighbour,
//            patchTypes,
//            patchNames,
//            patchSizes,
//            patchStarts,
//            procNeighbours,
//            cellCentres
//        );
    }
}


Foam::autoPtr<Foam::fvMesh> Foam::conformalVoronoiMesh::createDummyMesh
(
    const IOobject& io,
    const wordList& patchTypes,
    const wordList& patchNames,
    const labelList& patchSizes,
    const labelList& patchStarts,
    const labelList& procNeighbours
) const
{
    autoPtr<fvMesh> meshPtr
    (
        new fvMesh
        (
            io,
            xferCopy(pointField()),
            xferCopy(faceList()),
            xferCopy(cellList())
        )
    );
    fvMesh& mesh = meshPtr();

    List<polyPatch*> patches(patchStarts.size());

    forAll(patches, patchI)
    {
        if (patchTypes[patchI] == processorPolyPatch::typeName)
        {
            patches[patchI] = new processorPolyPatch
            (
                patchNames[patchI],
                0,          //patchSizes[p],
                0,          //patchStarts[p],
                patchI,
                mesh.boundaryMesh(),
                Pstream::myProcNo(),
                procNeighbours[patchI],
                coupledPolyPatch::COINCIDENTFULLMATCH
            );
        }
        else
        {
            patches[patchI] = polyPatch::New
            (
                patchTypes[patchI],
                patchNames[patchI],
                0,          //patchSizes[p],
                0,          //patchStarts[p],
                patchI,
                mesh.boundaryMesh()
            ).ptr();
        }
    }
    mesh.addFvPatches(patches);

    return meshPtr;
}


void Foam::conformalVoronoiMesh::checkProcessorPatchesMatch
(
    const wordList& patchTypes,
    const labelList& patchSizes,
    const labelList& procNeighbours
) const
{
    // Check patch sizes
    labelListList procPatchSizes
    (
        Pstream::nProcs(),
        labelList(Pstream::nProcs(), -1)
    );

    forAll(patchTypes, patchI)
    {
        if (patchTypes[patchI] == processorPolyPatch::typeName)
        {
            procPatchSizes[Pstream::myProcNo()][procNeighbours[patchI]]
                = patchSizes[patchI];
        }
    }

    Pstream::gatherList(procPatchSizes);

    if (Pstream::master())
    {
        bool allMatch = true;

        forAll(procPatchSizes, procI)
        {
            const labelList& patchSizes = procPatchSizes[procI];

            forAll(patchSizes, patchI)
            {
                if (patchSizes[patchI] != procPatchSizes[patchI][procI])
                {
                    allMatch = false;

                    Info<< indent << "Patches " << procI << " and " << patchI
                        << " have different sizes: " << patchSizes[patchI]
                        << " and " << procPatchSizes[patchI][procI] << endl;
                }
            }
        }

        if (allMatch)
        {
            Info<< indent << "All processor patches have matching numbers of "
                << "faces" << endl;
        }
    }
}


void Foam::conformalVoronoiMesh::reorderPoints
(
    pointField& points,
    labelList& boundaryPts,
    faceList& faces,
    const label nInternalFaces
) const
{
    Info<< incrIndent << indent << "Reordering points into internal/external"
        << endl;

    labelList oldToNew(points.size(), 0);

    // Find points that are internal
    for (label fI = nInternalFaces; fI < faces.size(); ++fI)
    {
        const face& f = faces[fI];

        forAll(f, fpI)
        {
            oldToNew[f[fpI]] = 1;
        }
    }

    const label nInternalPoints = points.size() - sum(oldToNew);

    label countInternal = 0;
    label countExternal = nInternalPoints;

    forAll(points, pI)
    {
        if (oldToNew[pI] == 0)
        {
            oldToNew[pI] = countInternal++;
        }
        else
        {
            oldToNew[pI] = countExternal++;
        }
    }

    Info<< indent
        << "Number of internal points: " << countInternal << nl
        << indent << "Number of external points: " << countExternal
        << decrIndent << endl;

    inplaceReorder(oldToNew, points);
    inplaceReorder(oldToNew, boundaryPts);

    forAll(faces, fI)
    {
        face& f = faces[fI];

        forAll(f, fpI)
        {
            f[fpI] = oldToNew[f[fpI]];
        }
    }
}


void Foam::conformalVoronoiMesh::reorderProcessorPatches
(
    const word& meshName,
    const fileName& instance,
    const pointField& points,
    faceList& faces,
    const wordList& patchTypes,
    const wordList& patchNames,
    const labelList& patchSizes,
    const labelList& patchStarts,
    const labelList& procNeighbours
) const
{
    Info<< incrIndent << indent << "Reordering processor patches" << endl;

    Info<< incrIndent;
    checkProcessorPatchesMatch(patchTypes, patchSizes, procNeighbours);

    // Create dummy mesh with correct proc boundaries to do sorting
    autoPtr<fvMesh> sortMeshPtr
    (
        createDummyMesh
        (
            IOobject
            (
                meshName,
                instance,
                runTime_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            patchTypes,
            patchNames,
            patchSizes,
            patchStarts,
            procNeighbours
        )
    );
    const fvMesh& sortMesh = sortMeshPtr();

    // Change the transform type on processors to coincident full match.
//    forAll(sortMesh.boundaryMesh(), patchI)
//    {
//        const polyPatch& patch = sortMesh.boundaryMesh()[patchI];
//
//        if (isA<processorPolyPatch>(patch))
//        {
//            const processorPolyPatch& cpPatch
//                = refCast<const processorPolyPatch>(patch);
//
//            processorPolyPatch& pPatch
//                = const_cast<processorPolyPatch&>(cpPatch);
//
//            pPatch.transform() = coupledPolyPatch::COINCIDENTFULLMATCH;
//        }
//    }

    // Rotation on new faces.
    labelList rotation(faces.size(), -1);
    labelList faceMap(faces.size(), -1);

    PstreamBuffers pBufs(Pstream::nonBlocking);

    // Send ordering
    forAll(sortMesh.boundaryMesh(), patchI)
    {
        const polyPatch& pp = sortMesh.boundaryMesh()[patchI];

        if (isA<processorPolyPatch>(pp))
        {
            refCast<const processorPolyPatch>(pp).initOrder
            (
                pBufs,
                primitivePatch
                (
                    SubList<face>
                    (
                        faces,
                        patchSizes[patchI],
                        patchStarts[patchI]
                    ),
                    points
                )
            );
        }
    }

    pBufs.finishedSends();

    // Receive and calculate ordering
    bool anyChanged = false;

    forAll(sortMesh.boundaryMesh(), patchI)
    {
        const polyPatch& pp = sortMesh.boundaryMesh()[patchI];

        if (isA<processorPolyPatch>(pp))
        {
            labelList patchFaceMap(patchSizes[patchI], -1);
            labelList patchFaceRotation(patchSizes[patchI], 0);

            bool changed = refCast<const processorPolyPatch>(pp).order
            (
                pBufs,
                primitivePatch
                (
                    SubList<face>
                    (
                        faces,
                        patchSizes[patchI],
                        patchStarts[patchI]
                    ),
                    points
                ),
                patchFaceMap,
                patchFaceRotation
            );

            if (changed)
            {
                // Merge patch face reordering into mesh face reordering table
                label start = patchStarts[patchI];

                forAll(patchFaceRotation, patchFaceI)
                {
                    rotation[patchFaceI + start]
                        = patchFaceRotation[patchFaceI];
                }

                forAll(patchFaceMap, patchFaceI)
                {
                    if (patchFaceMap[patchFaceI] != patchFaceI)
                    {
                        faceMap[patchFaceI + start]
                            = patchFaceMap[patchFaceI] + start;
                    }
                }

                anyChanged = true;
            }
        }
    }

    reduce(anyChanged, orOp<bool>());

    if (anyChanged)
    {
        inplaceReorder(faceMap, faces);

        // Rotate faces (rotation is already in new face indices).
        label nRotated = 0;

        forAll(rotation, faceI)
        {
            if (rotation[faceI] == -1)
            {
                continue;
            }

            if (rotation[faceI] != 0)
            {
                inplaceRotateList<List, label>(faces[faceI], rotation[faceI]);
                nRotated++;
            }
        }

        Info<< indent << returnReduce(nRotated, sumOp<label>())
            << " faces have been rotated" << decrIndent << decrIndent << endl;
    }
}


void Foam::conformalVoronoiMesh::writeMesh
(
    const word& meshName,
    const fileName& instance,
    pointField& points,
    labelList& boundaryPts,
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    const wordList& patchTypes,
    const wordList& patchNames,
    const labelList& patchSizes,
    const labelList& patchStarts,
    const labelList& procNeighbours,
    const pointField& cellCentres,
    const PackedBoolList& boundaryFacesToRemove
) const
{
    if (cvMeshControls().objOutput())
    {
        writeObjMesh(points, faces, word(meshName + ".obj"));
    }

    reorderPoints(points, boundaryPts, faces, patchStarts[0]);

    if (Pstream::parRun())
    {
        reorderProcessorPatches
        (
            meshName,
            instance,
            points,
            faces,
            patchTypes,
            patchNames,
            patchSizes,
            patchStarts,
            procNeighbours
        );
    }

    Info<< "    Constructing mesh" << endl;

    timeCheck("Before fvMesh construction");

    fvMesh mesh
    (
        IOobject
        (
            meshName,
            instance,
            runTime_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        xferMove(points),
        xferMove(faces),
        xferMove(owner),
        xferMove(neighbour)
    );

    Info<< "    Adding patches to mesh" << endl;

    List<polyPatch*> patches(patchStarts.size());

    label nValidPatches = 0;

    forAll(patches, p)
    {
        if (patchTypes[p] == processorPolyPatch::typeName)
        {
            // Do not create empty processor patches
            if (patchSizes[p] > 0)
            {
                patches[nValidPatches] = new processorPolyPatch
                (
                    patchNames[p],
                    patchSizes[p],
                    patchStarts[p],
                    nValidPatches,
                    mesh.boundaryMesh(),
                    Pstream::myProcNo(),
                    procNeighbours[p],
                    coupledPolyPatch::NOORDERING
                );

                nValidPatches++;
            }
        }
        else
        {
            patches[nValidPatches] = polyPatch::New
            (
                patchTypes[p],
                patchNames[p],
                patchSizes[p],
                patchStarts[p],
                nValidPatches,
                mesh.boundaryMesh()
            ).ptr();

            nValidPatches++;
        }
    }

    // Add indirectPatchFaces to a face zone
    {
        labelList addr(boundaryFacesToRemove.count());
        label count = 0;

        forAll(boundaryFacesToRemove, faceI)
        {
            if (boundaryFacesToRemove[faceI])
            {
                addr[count++] = faceI;
            }
        }

        label sz = mesh.faceZones().size();
        boolList flip(addr.size(), false);
        mesh.faceZones().setSize(sz + 1);
        mesh.faceZones().set
        (
            sz,
            new faceZone
            (
                "indirectPatchFaces",
                addr,
                flip,
                sz,
                mesh.faceZones()
            )
        );
    }

    patches.setSize(nValidPatches);

    mesh.addFvPatches(patches);

    timeCheck("Before fvMesh filtering");

    autoPtr<polyMeshFilter> meshFilter;

    label nInitialBadFaces = 0;

    if (cvMeshControls().filterEdges())
    {
        Info<< nl << "Filtering edges on polyMesh" << nl << endl;

        meshFilter.reset(new polyMeshFilter(mesh));

        // Filter small edges only. This reduces the number of faces so that
        // the face filtering is sped up.
        nInitialBadFaces = meshFilter().filterEdges(0);
        {
            const autoPtr<fvMesh>& newMesh = meshFilter().filteredMesh();

            polyTopoChange meshMod(newMesh);

            meshMod.changeMesh(mesh, false);
        }
    }

    if (cvMeshControls().filterFaces())
    {
        Info<< nl << "Filtering faces on polyMesh" << nl << endl;

        meshFilter.reset(new polyMeshFilter(mesh));

        meshFilter().filter(nInitialBadFaces);
        {
            const autoPtr<fvMesh>& newMesh = meshFilter().filteredMesh();

            polyTopoChange meshMod(newMesh);

            meshMod.changeMesh(mesh, false);
        }
    }

    timeCheck("After fvMesh filtering");

    mesh.setInstance(instance);

    if (!mesh.write())
    {
        FatalErrorIn("Foam::conformalVoronoiMesh::writeMesh(..)")
            << "Failed writing polyMesh."
            << exit(FatalError);
    }
    else
    {
        Info<< nl << "Written filtered mesh to "
            << mesh.polyMesh::instance() << nl
            << endl;
    }


    volTensorField alignments
    (
        IOobject
        (
            "alignmentsField",
            runTime_.timeName(),
            runTime_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        tensor::zero
    );

    forAll(mesh.cellCentres(), pI)
    {
        Vertex_handle nearV =
            nearest_vertex
            (
                toPoint<Point>(mesh.cellCentres()[pI])
            );
        alignments[pI] = nearV->alignment();
    }
    alignments.write();

    {
        volVectorField alignmentx
        (
            IOobject
            (
                "alignmentsx",
                runTime_.timeName(),
                runTime_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            vector::zero
        );
        forAll(alignmentx, aI)
        {
            alignmentx[aI] = alignments[aI].x();
        }
        alignmentx.write();
    }
    {
        volVectorField alignmenty
        (
            IOobject
            (
                "alignmentsy",
                runTime_.timeName(),
                runTime_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            vector::zero
        );
        forAll(alignmenty, aI)
        {
            alignmenty[aI] = alignments[aI].y();
        }
        alignmenty.write();
    }
    {
        volVectorField alignmentz
        (
            IOobject
            (
                "alignmentsz",
                runTime_.timeName(),
                runTime_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            vector::zero
        );
        forAll(alignmentz, aI)
        {
            alignmentz[aI] = alignments[aI].z();
        }
        alignmentz.write();
    }
    labelIOList boundaryIOPts
    (
        IOobject
        (
            "boundaryPoints",
            instance,
            runTime_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        boundaryPts
    );




    // Dump list of boundary points
    forAll(mesh.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        if (!isA<coupledPolyPatch>(pp))
        {
            forAll(pp, fI)
            {
                const face& boundaryFace = pp[fI];

                forAll(boundaryFace, pI)
                {
                    const label boundaryPointI = boundaryFace[pI];

                    boundaryIOPts[boundaryPointI] = boundaryPts[boundaryPointI];
                }
            }
        }
    }

    boundaryIOPts.write();

//    forAllConstIter(labelHashSet, pointsInPatch, pI)
//    {
//        const Foam::point& ptMaster = mesh.points()[pI.key()];
//
//        forAllConstIter(labelHashSet, pointsInPatch, ptI)
//        {
//            if (ptI.key() != pI.key())
//            {
//                const Foam::point& ptSlave = mesh.points()[ptI.key()];
//
//                const scalar dist = mag(ptMaster - ptSlave);
//                if (ptMaster == ptSlave)
//                {
//                    Pout<< "Point(" << pI.key() << ") " << ptMaster
//                        << " == "
//                        << "(" << ptI.key() << ") " << ptSlave
//                        << endl;
//                }
//                else if (dist == 0)
//                {
//                    Pout<< "Point(" << pI.key() << ") " << ptMaster
//                        << " ~= "
//                        << "(" << ptI.key() << ") " << ptSlave
//                        << endl;
//                }
//            }
//        }
//    }

//    writeCellSizes(mesh);

//    writeCellAlignments(mesh);

//    writeCellCentres(mesh);

//    findRemainingProtrusionSet(mesh);
}


void Foam::conformalVoronoiMesh::writeObjMesh
(
    const pointField& points,
    const faceList& faces,
    const fileName& fName
) const
{
    OFstream str(runTime_.path()/fName);

    Pout<< nl << "Writing points and faces to " << str.name() << endl;

    forAll(points, p)
    {
        meshTools::writeOBJ(str, points[p]);
    }

    forAll(faces, f)
    {
        str<< 'f';

        const face& fP = faces[f];

        forAll(fP, p)
        {
            str<< ' ' << fP[p] + 1;
        }

        str<< nl;
    }
}


void Foam::conformalVoronoiMesh::writeCellSizes
(
    const fvMesh& mesh
) const
{
    {
        timeCheck("Start writeCellSizes");

        Info<< nl << "Create targetCellSize volScalarField" << endl;

        volScalarField targetCellSize
        (
            IOobject
            (
                "targetCellSize",
                mesh.polyMesh::instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("cellSize", dimLength, 0),
            zeroGradientFvPatchScalarField::typeName
        );

        scalarField& cellSize = targetCellSize.internalField();

        const vectorField& C = mesh.cellCentres();

        forAll(cellSize, i)
        {
            cellSize[i] = cellShapeControls().cellSize(C[i]);
        }

        // Info<< nl << "Create targetCellVolume volScalarField" << endl;

        // volScalarField targetCellVolume
        // (
        //     IOobject
        //     (
        //         "targetCellVolume",
        //         mesh.polyMesh::instance(),
        //         mesh,
        //         IOobject::NO_READ,
        //         IOobject::AUTO_WRITE
        //     ),
        //     mesh,
        //     dimensionedScalar("cellVolume", dimLength, 0),
        //     zeroGradientFvPatchScalarField::typeName
        // );

        // targetCellVolume.internalField() = pow3(cellSize);

        // Info<< nl << "Create actualCellVolume volScalarField" << endl;

        // volScalarField actualCellVolume
        // (
        //     IOobject
        //     (
        //         "actualCellVolume",
        //         mesh.polyMesh::instance(),
        //         mesh,
        //         IOobject::NO_READ,
        //         IOobject::AUTO_WRITE
        //     ),
        //     mesh,
        //     dimensionedScalar("cellVolume", dimVolume, 0),
        //     zeroGradientFvPatchScalarField::typeName
        // );

        // actualCellVolume.internalField() = mesh.cellVolumes();

        // Info<< nl << "Create equivalentCellSize volScalarField" << endl;

        // volScalarField equivalentCellSize
        // (
        //     IOobject
        //     (
        //         "equivalentCellSize",
        //         mesh.polyMesh::instance(),
        //         mesh,
        //         IOobject::NO_READ,
        //         IOobject::AUTO_WRITE
        //     ),
        //     mesh,
        //     dimensionedScalar("cellSize", dimLength, 0),
        //     zeroGradientFvPatchScalarField::typeName
        // );

        // equivalentCellSize.internalField() = pow
        // (
        //     actualCellVolume.internalField(),
        //     1.0/3.0
        // );

        targetCellSize.correctBoundaryConditions();
        // targetCellVolume.correctBoundaryConditions();
        // actualCellVolume.correctBoundaryConditions();
        // equivalentCellSize.correctBoundaryConditions();

        targetCellSize.write();
        // targetCellVolume.write();
        // actualCellVolume.write();
        // equivalentCellSize.write();
    }

    // {
    //     polyMesh tetMesh
    //     (
    //         IOobject
    //         (
    //             "tetDualMesh",
    //             runTime_.constant(),
    //             runTime_,
    //             IOobject::MUST_READ
    //         )
    //     );

    //     pointMesh ptMesh(tetMesh);

    //     pointScalarField ptTargetCellSize
    //     (
    //         IOobject
    //         (
    //             "ptTargetCellSize",
    //             runTime_.timeName(),
    //             tetMesh,
    //             IOobject::NO_READ,
    //             IOobject::AUTO_WRITE
    //         ),
    //         ptMesh,
    //         dimensionedScalar("ptTargetCellSize", dimLength, 0),
    //         pointPatchVectorField::calculatedType()
    //     );

    //     scalarField& cellSize = ptTargetCellSize.internalField();

    //     const vectorField& P = tetMesh.points();

    //     forAll(cellSize, i)
    //     {
    //         cellSize[i] = cellShapeControls().cellSize(P[i]);
    //     }

    //     ptTargetCellSize.write();
    // }
}


void Foam::conformalVoronoiMesh::writeCellAlignments
(
    const fvMesh& mesh
) const
{
//    Info<< nl << "Create cellAlignments volTensorField" << endl;
//
//    volTensorField cellAlignments
//    (
//        IOobject
//        (
//            "cellAlignments",
//            mesh.polyMesh::instance(),
//            mesh,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        mesh,
//        tensor::I,
//        zeroGradientFvPatchTensorField::typeName
//    );
//
//    tensorField& cellAlignment = cellAlignments.internalField();
//
//    const vectorField& C = mesh.cellCentres();
//
//    vectorField xDir(cellAlignment.size());
//    vectorField yDir(cellAlignment.size());
//    vectorField zDir(cellAlignment.size());
//
//    forAll(cellAlignment, i)
//    {
//        cellAlignment[i] = cellShapeControls().cellAlignment(C[i]);
//        xDir[i] = cellAlignment[i] & vector(1, 0, 0);
//        yDir[i] = cellAlignment[i] & vector(0, 1, 0);
//        zDir[i] = cellAlignment[i] & vector(0, 0, 1);
//    }
//
//    OFstream xStr("xDir.obj");
//    OFstream yStr("yDir.obj");
//    OFstream zStr("zDir.obj");
//
//    forAll(xDir, i)
//    {
//        meshTools::writeOBJ(xStr, C[i], C[i] + xDir[i]);
//        meshTools::writeOBJ(yStr, C[i], C[i] + yDir[i]);
//        meshTools::writeOBJ(zStr, C[i], C[i] + zDir[i]);
//    }
//
//    cellAlignments.correctBoundaryConditions();
//
//    cellAlignments.write();
}


void Foam::conformalVoronoiMesh::writeCellCentres
(
    const fvMesh& mesh
) const
{
    Info<< "Writing components of cellCentre positions to volScalarFields"
        << " ccx, ccy, ccz in " <<  runTime_.timeName() << endl;

    for (direction i=0; i<vector::nComponents; i++)
    {
        volScalarField cci
        (
            IOobject
            (
                "cc" + word(vector::componentNames[i]),
                runTime_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh.C().component(i)
        );

        cci.write();
    }
}


Foam::labelHashSet Foam::conformalVoronoiMesh::findRemainingProtrusionSet
(
    const polyMesh& mesh
) const
{
    timeCheck("Start findRemainingProtrusionSet");

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    labelHashSet protrudingBoundaryPoints;

    forAll(patches, patchI)
    {
        const polyPatch& patch = patches[patchI];

        forAll(patch.localPoints(), pLPI)
        {
            label meshPtI = patch.meshPoints()[pLPI];

            const Foam::point& pt = patch.localPoints()[pLPI];

            if
            (
                geometryToConformTo_.wellOutside
                (
                    pt,
                    sqr(targetCellSize(pt))
                )
            )
            {
                protrudingBoundaryPoints.insert(meshPtI);
            }
        }
    }

    cellSet protrudingCells
    (
        mesh,
        "cvMesh_remainingProtrusions",
        mesh.nCells()/1000
    );

    forAllConstIter(labelHashSet, protrudingBoundaryPoints, iter)
    {
        const label pointI = iter.key();
        const labelList& pCells = mesh.pointCells()[pointI];

        forAll(pCells, pCI)
        {
            protrudingCells.insert(pCells[pCI]);
        }
    }

    label protrudingCellsSize = protrudingCells.size();

    reduce(protrudingCellsSize, sumOp<label>());

    if (cvMeshControls().objOutput() && protrudingCellsSize > 0)
    {
        Info<< nl << "Found " << protrudingCellsSize
            << " cells protruding from the surface, writing cellSet "
            << protrudingCells.name()
            << endl;

        protrudingCells.write();
    }

    return protrudingCells;
}


// ************************************************************************* //
