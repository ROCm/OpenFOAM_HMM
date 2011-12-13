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
#include "IOstreams.H"
#include "OFstream.H"
#include "pointMesh.H"
#include "pointFields.H"

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
            label mSize = m.size();
            label mPeak = m.peak();
            label mRss = m.rss();

            if (Pstream::parRun())
            {
                labelList allMSize(Pstream::nProcs());
                labelList allMPeak(Pstream::nProcs());
                labelList allMRss(Pstream::nProcs());

                allMSize[Pstream::myProcNo()] = mSize;
                allMPeak[Pstream::myProcNo()] = mPeak;
                allMRss[Pstream::myProcNo()] = mRss;

                Pstream::gatherList(allMSize);
                Pstream::gatherList(allMPeak);
                Pstream::gatherList(allMRss);

                Info<< "--- [ "
                    << "mem (kB) " << tab
                    << "size" << tab
                    << "peak" << tab
                    << "rss"
                    << " ] --- " << endl;

                forAll(allMSize, procI)
                {
                    Info<< "--- [ "
                        << procI << " " << tab
                        << allMSize[procI] << tab
                        << allMPeak[procI] << tab
                        << allMRss[procI]
                        << " ] --- " << endl;
                }

                Info<< "--- [ "
                    << "sum " << tab
                    << sum(allMSize) << tab
                    << sum(allMPeak) << tab
                    << sum(allMRss)
                    << " ] --- " << endl;

            }
            else
            {
                Info<< "--- [ "
                    << "mem size " << mSize << " kB, "
                    << "mem peak " << mPeak << " kB, "
                    << "mem rss " << mRss << " kB"
                    << " ] --- " << endl;
            }
        }
    }
}

void Foam::conformalVoronoiMesh::drawDelaunayCell
(
    Ostream& os,
    const Cell_handle& c,
    label offset
) const
{
    // Supply offset as tet number
    offset *= 5;

    os  << "# cell index: " << label(c->cellIndex()) << endl;

    os  << "# circumradius "
        << mag(topoint(dual(c)) - topoint(c->vertex(0)->point()))
        << endl;

    for (int i = 0; i < 4; i++)
    {
        os  << "# index type: "
            << label(c->vertex(i)->index()) << " "
            << label(c->vertex(i)->type()) << endl;

        meshTools::writeOBJ(os, topoint(c->vertex(i)->point()));
    }

    os  << "f " << 1 + offset << " " << 3 + offset << " " << 2 + offset << nl
        << "f " << 2 + offset << " " << 3 + offset << " " << 4 + offset << nl
        << "f " << 1 + offset << " " << 4 + offset << " " << 3 + offset << nl
        << "f " << 1 + offset << " " << 2 + offset << " " << 4 + offset << endl;

    os  << "# cicumcentre " << endl;

    meshTools::writeOBJ(os, topoint(dual(c)));

    os  << "l " << 1 + offset << " " << 5 + offset << endl;
}


void Foam::conformalVoronoiMesh::writePoints
(
    const fileName& fName,
    bool internalOnly
) const
{
    OFstream str(runTime_.path()/fName);

    Pout<< nl << "Writing points to " << str.name() << endl;

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (!internalOnly || vit->internalOrBoundaryPoint())
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


void Foam::conformalVoronoiMesh::writeMesh
(
    const fileName& instance,
    bool filterFaces
)
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
        faceList faces;
        labelList owner;
        labelList neighbour;
        wordList patchTypes;
        wordList patchNames;
        labelList patchSizes;
        labelList procNeighbours;
        pointField cellCentres;

        calcDualMesh
        (
            points,
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
            filterFaces
        );

        Info<< nl << "Writing polyMesh to " << instance << endl;

        writeMesh
        (
            Foam::polyMesh::defaultRegion,
            instance,
            points,
            faces,
            owner,
            neighbour,
            patchTypes,
            patchNames,
            patchSizes,
            dualPatchStarts,
            procNeighbours,
            cellCentres
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

        Info<< nl << "Writing tetDualMesh to " << instance << endl;

        writeMesh
        (
            "tetDualMesh",
            instance,
            points,
            faces,
            owner,
            neighbour,
            patchTypes,
            patchNames,
            patchSizes,
            patchStarts,
            procNeighbours,
            cellCentres
        );
    }
}


void Foam::conformalVoronoiMesh::writeMesh
(
    const word& meshName,
    const fileName& instance,
    pointField& points,
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    const wordList& patchTypes,
    const wordList& patchNames,
    const labelList& patchSizes,
    const labelList& patchStarts,
    const labelList& procNeighbours,
    const pointField& cellCentres
) const
{
    if (cvMeshControls().objOutput())
    {
        writeObjMesh(points, faces, word(meshName + ".obj"));
    }

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
                    procNeighbours[p]
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

    patches.setSize(nValidPatches);

    mesh.addFvPatches(patches);

    if (!mesh.write())
    {
        FatalErrorIn("Foam::conformalVoronoiMesh::writeMesh(..)")
            << "Failed writing polyMesh."
            << exit(FatalError);
    }

    // pointIOField cellCs
    // (
    //     IOobject
    //     (
    //         "cellCentres",
    //         mesh.pointsInstance(),
    //         polyMesh::meshSubDir,
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     cellCentres
    // );

    // Info<< nl
    //     << "Writing " << cellCs.name()
    //     << " to " << cellCs.instance()
    //     << endl;

    // cellCs.write();

    writeCellSizes(mesh);

    writeCellCentres(mesh);

    findRemainingProtrusionSet(mesh);
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
            cellSize[i] = cellSizeControl().cellSize(C[i]);
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
    //         cellSize[i] = cellSizeControl().cellSize(P[i]);
    //     }

    //     ptTargetCellSize.write();
    // }
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


void Foam::conformalVoronoiMesh::findRemainingProtrusionSet
(
    const fvMesh& mesh
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

    if (protrudingCellsSize > 0)
    {
        Info<< nl << "Found " << protrudingCellsSize
            << " cells protruding from the surface, writing cellSet "
            << protrudingCells.name()
            << endl;

        protrudingCells.write();
    }
}


// ************************************************************************* //
