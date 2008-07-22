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

\*----------------------------------------------------------------------------*/

#include "CV3D.H"
#include "plane.H"
#include "triSurfaceTools.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::CV3D::insertFeaturePoints()
{
    labelList featPoints(qSurf_.extractFeatures3D(controls_.featAngle));

    const pointField& localPts = qSurf_.localPoints();

    forAll(featPoints, i)
    {
        label ptI = featPoints[i];

        const point& featPt = localPts[ptI];

        insertPoint(featPt, Vb::MIRROR_POINT);
    }

    if (controls_.writeFeatureTriangulation)
    {
        writePoints("feat_allPoints.obj", false);
        writePoints("feat_points.obj", true);
//         writeFaces("feat_allFaces.obj", false);
//         writeFaces("feat_faces.obj", true);
        writeTriangles("feat_triangles.obj", true);
    }
}


// ************************************************************************* //
