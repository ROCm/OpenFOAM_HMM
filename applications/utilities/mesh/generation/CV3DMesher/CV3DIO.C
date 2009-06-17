/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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

#include "CV3D.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CV3D::writePoints(const fileName& fName, bool internalOnly) const
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


void Foam::CV3D::writeDual
(
    const pointField& points,
    const faceList& faces,
    const fileName& fName
) const
{
    Info<< nl << "Writing dual points and faces to " << fName << endl;

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


void Foam::CV3D::writeTriangles(const fileName& fName, bool internalOnly) const
{
    Info<< nl << "Writing triangles to " << fName << endl;
    OFstream str(fName);

    labelList vertexMap(number_of_vertices());
    label verti = 0;

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (!internalOnly || !vit->farPoint())
        {
            vertexMap[vit->index()] = verti++;
            meshTools::writeOBJ(str, topoint(vit->point()));
        }
    }

    for
    (
        Triangulation::Finite_facets_iterator fit = finite_facets_begin();
        fit != finite_facets_end();
        ++fit
    )
    {
        const Cell_handle& c(fit->first);

        const int& oppositeVertex(fit->second);

        List<label> facetIndices(3,-1);

        bool writeFacet = true;

        for(label i = 0, k = 0;i < 4; i++)
        {
            if(i != oppositeVertex)
            {
                if(!internalOnly || !c->vertex(i)->farPoint())
                {
                    facetIndices[k] = i;
                    k++;
                }
                else
                {
                    writeFacet = false;
                }
            }
        }

        if(writeFacet)
        {
            str << "f "
                << vertexMap[c->vertex(facetIndices[0])->index()] + 1
                << ' ' << vertexMap[c->vertex(facetIndices[1])->index()] + 1
                << ' ' << vertexMap[c->vertex(facetIndices[2])->index()] + 1
                << nl;
        }
    }
}

void Foam::CV3D::writeMesh(bool writeToTimestep)
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

    writeDual(points, faces, "dualMesh.obj");

    IOobject io
    (
        Foam::polyMesh::defaultRegion,
        runTime_.constant(),
        runTime_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    );

    if (writeToTimestep)
    {
        Info<< nl << "Writing polyMesh to time directory "
            << runTime_.timeName() << endl;

        io = IOobject
        (
            Foam::polyMesh::defaultRegion,
            runTime_.path()/runTime_.timeName(),
            runTime_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        );
    }
    else
    {
        Info<< nl << "Writing polyMesh to constant." << endl;
    }

    polyMesh pMesh
    (
        io,
        xferMove(points),
        xferMove(faces),
        xferMove(owner),
        xferMove(neighbour)
    );
    // polyMesh pMesh
    // (
    //     io,
    //     points,
    //     faces,
    //     owner,
    //     neighbour
    // );

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
        FatalErrorIn("CV3D::writeMesh()")
            << "Failed writing polyMesh."
            << exit(FatalError);
    }


}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// ************************************************************************* //
