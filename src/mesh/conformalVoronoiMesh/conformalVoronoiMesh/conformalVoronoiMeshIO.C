/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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
#include "IOstreams.H"
#include "OFstream.H"
#include "zeroGradientPointPatchField.H"
#include "pointMesh.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::conformalVoronoiMesh::writePoints
(
    const fileName& fName,
    bool internalOnly
) const
{
    Info<< nl << "Writing points to " << fName << endl;

    OFstream str(fName);

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
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
    const List<point>& points
) const
{
    if (points.size())
    {
        Info<< nl << "Writing " << points.size() << " points from pointList to "
            << fName << endl;

        OFstream str(fName);

        forAll(points, p)
        {
            meshTools::writeOBJ(str, points[p]);
        }
    }
}


void Foam::conformalVoronoiMesh::writeInternalDelaunayVertices
(
    bool writeToConstant
) const
{
    pointField internalDelaunayVertices(number_of_vertices());

    label vertI = 0;

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
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

    IOobject io
    (
        "internalDelaunayVertices",
        runTime_.timeName(),
        runTime_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    );

    if (writeToConstant)
    {
        io = IOobject
        (
            "internalDelaunayVertices",
            runTime_.constant(),
            runTime_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        );
    }

    Info<< nl << "Writing " << io.name() << " to " << io.instance() << endl;

    pointIOField internalDVs(io, internalDelaunayVertices);

    internalDVs.write();
}


void Foam::conformalVoronoiMesh::writeMesh(bool writeToConstant)
{
    writeInternalDelaunayVertices(writeToConstant);

    // {
    //     pointField points;
    //     faceList faces;
    //     labelList owner;
    //     labelList neighbour;
    //     wordList patchNames;
    //     labelList patchSizes;
    //     labelList patchStarts;

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

    //     writeMesh
    //     (
    //         "tetDualMesh",
    //         writeToConstant,
    //         points,
    //         faces,
    //         owner,
    //         neighbour,
    //         patchNames,
    //         patchSizes,
    //         patchStarts
    //     );
    // }

    {
        pointField points;
        faceList faces;
        labelList owner;
        labelList neighbour;
        wordList patchNames;
        labelList patchSizes;
        labelList patchStarts;

        calcDualMesh
        (
            points,
            faces,
            owner,
            neighbour,
            patchNames,
            patchSizes,
            patchStarts
        );

        writeMesh
        (
            Foam::polyMesh::defaultRegion,
            writeToConstant,
            points,
            faces,
            owner,
            neighbour,
            patchNames,
            patchSizes,
            patchStarts
        );
    }
}

void Foam::conformalVoronoiMesh::writeMesh
(
    const word& meshName,
    bool writeToConstant,
    pointField& points,
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    wordList& patchNames,
    labelList& patchSizes,
    labelList& patchStarts
) const
{
    if(cvMeshControls().objOutput())
    {
        writeObjMesh(points, faces, word(meshName + ".obj"));
    }

    IOobject io
    (
        meshName,
        runTime_.timeName(),
        runTime_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    );

    if (writeToConstant)
    {
        Info<< nl << "Writing mesh to constant." << endl;

        io = IOobject
        (
            meshName,
            runTime_.constant(),
            runTime_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        );
    }
    else
    {
        Info<< nl << "Writing mesh to time directory "
            << runTime_.timeName() << endl;
    }

    polyMesh pMesh
    (
        io,
        xferMove(points),
        xferMove(faces),
        xferMove(owner),
        xferMove(neighbour)
    );

    List<polyPatch*> patches(patchStarts.size());

    forAll (patches, p)
    {
        patches[p] = new polyPatch
        (
            patchNames[p],
            patchSizes[p],
            patchStarts[p],
            p,
            pMesh.boundaryMesh()
        );
    }

    pMesh.addPatches(patches);

    if (!pMesh.write())
    {
        FatalErrorIn("Foam::conformalVoronoiMesh::writeMesh")
            << "Failed writing polyMesh."
            << exit(FatalError);
    }
}


void Foam::conformalVoronoiMesh::writeObjMesh
(
    const pointField& points,
    const faceList& faces,
    const fileName& fName
) const
{
    Info<< nl << "Writing points and faces to " << fName << endl;

    OFstream str(fName);

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


void Foam::conformalVoronoiMesh::writeTargetCellSize() const
{
    {
        Info<< nl << "Create fvMesh" << endl;

        fvMesh fMesh
        (
            IOobject
            (
                Foam::polyMesh::defaultRegion,
                runTime_.constant(),
                runTime_,
                IOobject::MUST_READ
            )
        );

        timeCheck();

        Info<< nl << "Create targetCellSize volScalarField" << endl;

        volScalarField targetCellSize
        (
            IOobject
            (
                "targetCellSize",
                runTime_.timeName(),
                runTime_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fMesh,
            dimensionedScalar("cellSize", dimLength, 0),
            zeroGradientPointPatchField<scalar>::typeName
        );

        scalarField& cellSize = targetCellSize.internalField();

        const vectorField& C = fMesh.cellCentres();

        forAll(cellSize, i)
        {
            cellSize[i] = cellSizeControl().cellSize(C[i]);
        }

        targetCellSize.write();
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


// ************************************************************************* //
