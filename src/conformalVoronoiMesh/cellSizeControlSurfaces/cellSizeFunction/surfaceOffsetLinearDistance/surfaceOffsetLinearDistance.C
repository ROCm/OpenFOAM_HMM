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

#include "surfaceOffsetLinearDistance.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(surfaceOffsetLinearDistance, 0);
addToRunTimeSelectionTable
(
    cellSizeFunction,
    surfaceOffsetLinearDistance,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

surfaceOffsetLinearDistance::surfaceOffsetLinearDistance
(
    const dictionary& initialPointsDict,
    const conformalVoronoiMesh& cvMesh,
    const searchableSurface& surface
)
:
    cellSizeFunction(typeName, initialPointsDict, cvMesh, surface),
    surfaceCellSize_(readScalar(coeffsDict().lookup("surfaceCellSize"))),
    distanceCellSize_(readScalar(coeffsDict().lookup("distanceCellSize"))),
    surfaceOffset_(readScalar(coeffsDict().lookup("surfaceOffset"))),
    totalDistance_
    (
        readScalar(coeffsDict().lookup("distance")) + surfaceOffset_
    ),
    totalDistanceSqr_(sqr(totalDistance_)),
    gradient_
    (
        (distanceCellSize_ - surfaceCellSize_)/(totalDistance_ - surfaceOffset_)
    ),
    intercept_(surfaceCellSize_ - gradient_*surfaceOffset_)
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

scalar surfaceOffsetLinearDistance::sizeFunction(scalar d) const
{
    if (d <= surfaceOffset_)
    {
        return surfaceCellSize_;
    }

    return gradient_*d + intercept_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool surfaceOffsetLinearDistance::cellSize(const point& pt, scalar& size) const
{
    size = 0;

    List<pointIndexHit> hits;

    surface_.findNearest
    (
        pointField(1, pt),
        scalarField(1, totalDistanceSqr_),
        hits
    );

    const pointIndexHit& hitInfo = hits[0];

    if (hitInfo.hit())
    {
        if (sideMode_ == BOTHSIDES)
        {
            size = sizeFunction(mag(pt - hitInfo.hitPoint()));

            return true;
        }

        pointField ptF(1, pt);
        List<searchableSurface::volumeType> vTL;

        surface_.getVolumeType(ptF, vTL);

        bool functionApplied = false;

        if
        (
            sideMode_ == INSIDE
         && vTL[0] == searchableSurface::INSIDE
        )
        {
            size = sizeFunction(mag(pt - hitInfo.hitPoint()));

            functionApplied = true;
        }
        else if
        (
            sideMode_ == OUTSIDE
         && vTL[0] == searchableSurface::OUTSIDE
        )
        {
            size = sizeFunction(mag(pt - hitInfo.hitPoint()));

            functionApplied = true;
        }

        return functionApplied;
    }

    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
