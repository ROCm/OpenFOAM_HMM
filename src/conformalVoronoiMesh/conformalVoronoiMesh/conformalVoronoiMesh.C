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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conformalVoronoiMesh::conformalVoronoiMesh
(
    const Time& runTime,
    const IOdictionary& cvMeshDict
)
:
    HTriangulation(),
    runTime_(runTime),
    cvMeshControls_(*this, cvMeshDict),
    startOfInternalPoints_(0),
    startOfSurfacePointPairs_(0),
    initialPointsMethod_
    (
        initialPointsMethod::New
        (
            cvMeshDict.subDict("surfaceConformation").subDict("initialPoints"),
            *this
        )
    )
{
    insertInitialPoints();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conformalVoronoiMesh::~conformalVoronoiMesh()
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::conformalVoronoiMesh::insertInitialPoints()
{
    startOfInternalPoints_ = number_of_vertices();

    label nVert = startOfInternalPoints_;

    Info<< nl << "Inserting initial points" << nl
        << "    " << nVert << " existing points" << endl;

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
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
