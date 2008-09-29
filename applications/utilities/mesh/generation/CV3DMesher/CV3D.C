/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
#include "Random.H"
#include "IFstream.H"
#include "uint.H"
#include "ulong.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::CV3D::insertBoundingBox()
{
    Info<< "insertBoundingBox: creating bounding mesh" << nl << endl;
    scalar bigSpan = 10*tols_.span;
    insertPoint(point(-bigSpan, -bigSpan, -bigSpan), Vb::FAR_POINT);
    insertPoint(point(-bigSpan, -bigSpan,  bigSpan), Vb::FAR_POINT);
    insertPoint(point(-bigSpan,  bigSpan, -bigSpan), Vb::FAR_POINT);
    insertPoint(point(-bigSpan,  bigSpan,  bigSpan), Vb::FAR_POINT);
    insertPoint(point( bigSpan, -bigSpan, -bigSpan), Vb::FAR_POINT);
    insertPoint(point( bigSpan, -bigSpan,  bigSpan), Vb::FAR_POINT);
    insertPoint(point( bigSpan,  bigSpan, -bigSpan), Vb::FAR_POINT);
    insertPoint(point( bigSpan,  bigSpan , bigSpan), Vb::FAR_POINT);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CV3D::CV3D
(
    const dictionary& controlDict,
    const querySurface& qSurf
)
:
    HTriangulation(),
    qSurf_(qSurf),
    controls_(controlDict),
    tols_(controlDict, controls_.minCellSize, qSurf.bb()),
    startOfInternalPoints_(0),
    startOfSurfacePointPairs_(0)
{
    // insertBoundingBox();
    insertFeaturePoints();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::CV3D::~CV3D()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CV3D::insertPoints
(
    const pointField& points,
    const scalar nearness
)
{
    Info<< "insertInitialPoints(const pointField& points): ";

    startOfInternalPoints_ = number_of_vertices();
    label nVert = startOfInternalPoints_;

    // Add the points and index them
    forAll(points, i)
    {
        const point& p = points[i];

        if (qSurf_.wellInside(p, nearness))
        {
            insert(toPoint(p))->index() = nVert++;
        }
        else
        {
            Warning
                << "Rejecting point " << p << " outside surface" << endl;
        }
    }

    Info<< nVert << " vertices inserted" << endl;

    if (controls_.writeInitialTriangulation)
    {
        // Checking validity of triangulation
        assert(is_valid());

        writeTriangles("initial_triangles.obj", true);
//         writeFaces("initial_faces.obj", true);
    }
}


void Foam::CV3D::insertPoints(const fileName& pointFileName)
{
    IFstream pointsFile(pointFileName);

    if (pointsFile.good())
    {
        insertPoints(pointField(pointsFile), 0.5*controls_.minCellSize2);
    }
    else
    {
        FatalErrorIn("insertInitialPoints")
            << "Could not open pointsFile " << pointFileName
            << exit(FatalError);
    }
}


void Foam::CV3D::insertGrid()
{
    Info<< "insertInitialGrid: ";

    startOfInternalPoints_ = number_of_vertices();
    label nVert = startOfInternalPoints_;

    scalar x0 = qSurf_.bb().min().x();
    scalar xR = qSurf_.bb().max().x() - x0;
    int ni = int(xR/controls_.minCellSize) + 1;

    scalar y0 = qSurf_.bb().min().y();
    scalar yR = qSurf_.bb().max().y() - y0;
    int nj = int(yR/controls_.minCellSize) + 1;

    scalar z0 = qSurf_.bb().min().z();
    scalar zR = qSurf_.bb().max().z() - z0;
    int nk = int(zR/controls_.minCellSize) + 1;

    vector delta(xR/ni, yR/nj, zR/nk);

    delta *= pow((1.0/2.0),-(1.0/3.0));

    Random rndGen(1321);
    scalar pert = controls_.randomPerturbation*cmptMin(delta);

    for (int i=0; i<ni; i++)
    {
        for (int j=0; j<nj; j++)
        {
            for (int k=0; k<nk; k++)
            {
                point p1
                (
                    x0 + i*delta.x(),
                    y0 + j*delta.y(),
                    z0 + k*delta.z()
                );

                point p2 = p1 + 0.5*delta;

                if (controls_.randomiseInitialGrid)
                {
                    p1.x() += pert*(rndGen.scalar01() - 0.5);
                    p1.y() += pert*(rndGen.scalar01() - 0.5);
                    p1.z() += pert*(rndGen.scalar01() - 0.5);
                }

                if (qSurf_.wellInside(p1, 0.5*controls_.minCellSize2))
                {
                    insert(Point(p1.x(), p1.y(), p1.z()))->index() = nVert++;
                }

                if (controls_.randomiseInitialGrid)
                {
                    p2.x() += pert*(rndGen.scalar01() - 0.5);
                    p2.y() += pert*(rndGen.scalar01() - 0.5);
                    p2.z() += pert*(rndGen.scalar01() - 0.5);
                }

                if (qSurf_.wellInside(p2, 0.5*controls_.minCellSize2))
                {
                    insert(Point(p2.x(), p2.y(), p2.z()))->index() = nVert++;
                }
            }
        }
    }

    Info<< nVert << " vertices inserted" << nl << endl;

    if (controls_.writeInitialTriangulation)
    {
//         Checking validity of triangulation
        assert(is_valid());

        writePoints("initial_points.obj", true);
        writeTriangles("initial_triangles.obj", true);
//         writeFaces("initial_faces.obj", true);
    }
}


void Foam::CV3D::relaxPoints(const scalar relaxation)
{
    Info<< "Calculating new points: " << nl << endl;

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalPoint())
        {
            // movePoint(vit, newPoint);
        }
    }

}


void Foam::CV3D::insertSurfacePointPairs()
{
    startOfSurfacePointPairs_ = number_of_vertices();

    if (controls_.insertSurfaceNearestPointPairs)
    {
        insertSurfaceNearestPointPairs();
    }

    if (controls_.writeNearestTriangulation)
    {
        // writeFaces("near_allFaces.obj", false);
        // writeFaces("near_faces.obj", true);
        writeTriangles("near_triangles.obj", true);
    }

    if (controls_.insertSurfaceNearPointPairs)
    {
        insertSurfaceNearPointPairs();
    }

    startOfBoundaryConformPointPairs_ = number_of_vertices();
}

void Foam::CV3D::boundaryConform()
{
}


void Foam::CV3D::removeSurfacePointPairs()
{
    Info<< "Removing surface point pairs." << nl << endl;

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->index() >= startOfSurfacePointPairs_)
        {
            remove(vit);
        }
    }
}


void Foam::CV3D::write() const
{
    if (controls_.writeFinalTriangulation)
    {
        writePoints("allPoints.obj", false);
        writePoints("points.obj", true);
//         writeFaces("allFaces.obj", false);
//         writeFaces("faces.obj", true);
        writeTriangles("allTriangles.obj", false);
        writeTriangles("triangles.obj", true);
//         writeMesh();
    }
}

// ************************************************************************* //
