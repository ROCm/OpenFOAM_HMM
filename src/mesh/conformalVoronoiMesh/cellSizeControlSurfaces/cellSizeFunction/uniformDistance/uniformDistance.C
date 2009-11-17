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

#include "uniformDistance.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(uniformDistance, 0);
addToRunTimeSelectionTable(cellSizeFunction, uniformDistance, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

uniformDistance::uniformDistance
(
    const dictionary& initialPointsDict,
    const conformalVoronoiMesh& cvMesh,
    const searchableSurface& surface
)
:
    cellSizeFunction(typeName, initialPointsDict, cvMesh, surface),
    cellSize_(readScalar(coeffsDict().lookup("cellSize"))),
    distance_(readScalar(coeffsDict().lookup("distance"))),
    distanceSqr_(sqr(distance_))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool uniformDistance::cellSize
(
    const point& pt,
    scalar& size,
    bool isSurfacePoint
) const
{
    if (isSurfacePoint)
    {
        size = cellSize_;

        return true;
    }

    size = 0;

    List<pointIndexHit> hits;

    surface_.findNearest
    (
        pointField(1, pt),
        scalarField(1, distanceSqr_),
        hits
    );

    const pointIndexHit& hitInfo = hits[0];

    if (hitInfo.hit())
    {
        if (sideMode_ == rmBothsides)
        {
            size = cellSize_;

            return true;
        }

        // If the nearest point is essentially on the surface, do not do a
        // getVolumeType calculation, as it will be prone to error.
        if (mag(pt  - hitInfo.hitPoint()) < snapToSurfaceTol_)
        {
            size = cellSize_;

            return true;
        }

        pointField ptF(1, pt);
        List<searchableSurface::volumeType> vTL;

        surface_.getVolumeType(ptF, vTL);

        bool functionApplied = false;

        if
        (
            sideMode_ == smInside
         && vTL[0] == searchableSurface::INSIDE
        )
        {
            size = cellSize_;

            functionApplied = true;
        }
        else if
        (
            sideMode_ == smOutside
         && vTL[0] == searchableSurface::OUTSIDE
        )
        {
            size = cellSize_;

            functionApplied = true;
        }

        return functionApplied;
    }

    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
