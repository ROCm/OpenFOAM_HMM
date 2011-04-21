/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
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
#include "zeroGradientPointPatchField.H"
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
                reduce(mSize, sumOp<label>());
                reduce(mPeak, sumOp<label>());
                reduce(mRss, sumOp<label>());
            }

            Info<< "--- [ "
                << runTime_.elapsedCpuTime() << " s, "
                << "mem size " << mSize << " kB, "
                << "mem peak " << mPeak << " kB, "
                << "mem rss " << mRss << " kB"
                << " ] --- " << endl;
        }
    }
}

void Foam::conformalVoronoiMesh::drawDelaunayCell
(
    Ostream& os,
    const Delaunay::Finite_cells_iterator& cit
) const
{
    for (int i = 0; i < 4; i++)
    {
        os  << "# index type: "
            << label(cit->vertex(i)->index()) << " "
            << label(cit->vertex(i)->type()) << endl;

        meshTools::writeOBJ(os, topoint(cit->vertex(i)->point()));
    }

    os  << "f 1 3 2" << nl
        << "f 2 3 4" << nl
        << "f 1 4 3" << nl
        << "f 1 2 4" << endl;

    meshTools::writeOBJ(os, topoint(dual(cit)));

    os  << "l 1 5" << endl;
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

    // {
    //     pointField points;
    //     faceList faces;
    //     labelList owner;
    //     labelList neighbour;
    //     wordList patchNames;
    //     labelList patchSizes;
    //     labelList patchStarts;
    //     pointField cellCentres;

    //     calcTetMesh
    //     (
    //         points,
    //         faces,
    //         owner,
    //         neighbour,
    //         patchNames,
    //         patchSizes,
    //         patchStarts
    //     );

    //     Info<< nl << "Writing tetDualMesh to " << instance << endl;

    //     writeMesh
    //     (
    //         "tetDualMesh",
    //         instance,
    //         points,
    //         faces,
    //         owner,
    //         neighbour,
    //         patchNames,
    //         patchSizes,
    //         patchStarts,
    //         cellCentres
    //     );
    // }

    {
        pointField points;
        faceList faces;
        labelList owner;
        labelList neighbour;
        wordList patchTypes;
        wordList patchNames;
        labelList patchSizes;
        labelList patchStarts;
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
            patchStarts,
            procNeighbours,
            cellCentres,
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
    wordList& patchTypes,
    wordList& patchNames,
    labelList& patchSizes,
    labelList& patchStarts,
    labelList& procNeighbours,
    const pointField& cellCentres
) const
{
    if(cvMeshControls().objOutput())
    {
        writeObjMesh(points, faces, word(meshName + ".obj"));
    }

    polyMesh mesh
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

    forAll(patches, p)
    {
        if (patchTypes[p] == processorPolyPatch::typeName)
        {
            patches[p] = new processorPolyPatch
            (
                patchNames[p],
                patchSizes[p],
                patchStarts[p],
                p,
                mesh.boundaryMesh(),
                Pstream::myProcNo(),
                procNeighbours[p]
            );
        }
        else
        {
            patches[p] = polyPatch::New
            (
                patchTypes[p],
                patchNames[p],
                patchSizes[p],
                patchStarts[p],
                p,
                mesh.boundaryMesh()
            ).ptr();
        }
    }

    Info<< "addPatches(patches, false); FALSE REQUIRED TO AVOID MEMORY CRASH"
        << endl;

    // mesh.addPatches(patches);
    mesh.addPatches(patches, false);

    if (!mesh.write())
    {
        FatalErrorIn("Foam::conformalVoronoiMesh::writeMesh")
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


    Info<< "DISABLED WRITING OF CELL SIZE AND PROTRUSION SET" << endl;

    // writeCellSizes(mesh);

    // findRemainingProtrusionSet(mesh);
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

    forAll (faces, f)
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
            zeroGradientPointPatchField<scalar>::typeName
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
        //     zeroGradientPointPatchField<scalar>::typeName
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
        //     zeroGradientPointPatchField<scalar>::typeName
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
        //     zeroGradientPointPatchField<scalar>::typeName
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

    if (!protrudingCells.empty())
    {
        Pout<< nl << "Found " << protrudingCells.size()
            << " cells protruding from the surface, writing cellSet "
            << protrudingCells.name()
            << endl;

        protrudingCells.write();
    }
}


// ************************************************************************* //
