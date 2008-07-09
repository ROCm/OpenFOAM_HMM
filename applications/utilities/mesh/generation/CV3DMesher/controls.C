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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CV3D::controls::controls(const dictionary& controlDict)
:
    minCellSize(readScalar(controlDict.lookup("minCellSize"))),
    minCellSize2(Foam::sqr(minCellSize)),
    featAngle(readScalar(controlDict.lookup("featureAngle"))),
    maxQuadAngle(readScalar(controlDict.lookup("maxQuadAngle"))),
    insertSurfaceNearestPointPairs
    (
        controlDict.lookup("insertSurfaceNearestPointPairs")
    ),
    insertSurfaceNearPointPairs
    (
        controlDict.lookup("insertSurfaceNearPointPairs")
    ),
    writeInitialTriangulation(controlDict.lookup("writeInitialTriangulation")),
    writeFeatureTriangulation(controlDict.lookup("writeFeatureTriangulation")),
    writeNearestTriangulation(controlDict.lookup("writeNearestTriangulation")),
    writeInsertedPointPairs(controlDict.lookup("writeInsertedPointPairs")),
    writeFinalTriangulation(controlDict.lookup("writeFinalTriangulation")),
    randomiseInitialGrid(controlDict.lookup("randomiseInitialGrid")),
    randomPerturbation(readScalar(controlDict.lookup("randomPerturbation")))
{}


// ************************************************************************* //
