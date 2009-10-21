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

#include "uniform.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(uniform, 0);
addToRunTimeSelectionTable(cellSizeFunction, uniform, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

uniform::uniform
(
    const dictionary& initialPointsDict,
    const conformalVoronoiMesh& cvMesh,
    const searchableSurface& surface
)
:
    cellSizeFunction(typeName, initialPointsDict, cvMesh, surface),
    cellSize_(readScalar(coeffsDict().lookup("cellSize")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool uniform::cellSize
(
    const point& pt,
    scalar& size,
    bool isSurfacePoint
) const
{
    if (sideMode_ == BOTHSIDES || isSurfacePoint)
    {
        size = cellSize_;

        return true;
    }

    size = 0;

    List<pointIndexHit> hits;

    surface_.findNearest
    (
        pointField(1, pt),
        scalarField(1, sqr(snapToSurfaceTol_)),
        hits
    );

    const pointIndexHit& hitInfo = hits[0];

    // If the nearest point is essentially on the surface, do not do a
    // getVolumeType calculation, as it will be prone to error.
    if (hitInfo.hit())
    {
        size = cellSize_;

        return true;
    }

    pointField ptF(1, pt);
    List<searchableSurface::volumeType> vTL(1);

    surface_.getVolumeType(ptF, vTL);

    bool functionApplied = false;

    if
    (
        sideMode_ == INSIDE
     && vTL[0] == searchableSurface::INSIDE
    )
    {
        size = cellSize_;

        functionApplied = true;
    }
    else if
    (
        sideMode_ == OUTSIDE
     && vTL[0] == searchableSurface::OUTSIDE
    )
    {
        size = cellSize_;

        functionApplied = true;
    }

    return functionApplied;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
