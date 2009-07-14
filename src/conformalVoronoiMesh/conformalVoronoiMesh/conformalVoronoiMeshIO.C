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


void Foam::conformalVoronoiMesh::writeInternalDelaunayVertices() const
{
    Info<< nl << "    Writing internal Delaunay vertices to pointField "
        << "ADD NAME CHOICE ARGUMENT" << endl;

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

    pointIOField internalDVs
    (
        IOobject
        (
            "internalDelaunayVertices",
            runTime_.timeName(),
            runTime_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        internalDelaunayVertices
    );

    internalDVs.write();
}


void Foam::conformalVoronoiMesh::writeMesh(bool writeToConstant)
{
    pointField points(0);
    faceList faces(0);
    labelList owner(0);
    labelList neighbour(0);
    wordList patchNames(0);
    labelList patchSizes(0);
    labelList patchStarts(0);

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

    if(cvMeshControls().objOutput())
    {
        writeDual(points, faces, "dualMesh.obj");
    }

    IOobject io
    (
        Foam::polyMesh::defaultRegion,
        runTime_.timeName(),
        runTime_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    );

    if (writeToConstant)
    {
        Info<< nl << "    Writing polyMesh to constant." << endl;

        io = IOobject
        (
            Foam::polyMesh::defaultRegion,
            runTime_.constant(),
            runTime_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        );

    }
    else
    {
        Info<< nl << "    Writing polyMesh to time directory "
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
        FatalErrorIn("Foam::conformalVoronoiMesh::writeMesh()")
        << "Failed writing polyMesh."
            << exit(FatalError);
    }
}


void Foam::conformalVoronoiMesh::writeDual
(
    const pointField& points,
    const faceList& faces,
    const fileName& fName
) const
{
    Info<< nl << "    Writing dual points and faces to " << fName << endl;

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


// ************************************************************************* //
