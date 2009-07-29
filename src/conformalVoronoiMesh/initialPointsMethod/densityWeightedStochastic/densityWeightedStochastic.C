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

#include "densityWeightedStochastic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(densityWeightedStochastic, 0);
addToRunTimeSelectionTable
(
    initialPointsMethod,
    densityWeightedStochastic,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

densityWeightedStochastic::densityWeightedStochastic
(
    const dictionary& initialPointsDict,
    const conformalVoronoiMesh& cvMesh
)
:
    initialPointsMethod(typeName, initialPointsDict, cvMesh),
    totalVolume_(readScalar(detailsDict().lookup("totalVolume"))),
    maxDensity_
    (
        1.0/pow3(readScalar(detailsDict().lookup("minCellSize")))
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

std::vector<Vb::Point> densityWeightedStochastic::initialPoints() const
{
    const boundBox& bb = cvMesh_.geometryToConformTo().bounds();

    Random rndGen(5234986);

    std::vector<Vb::Point> initialPoints;

    scalar volumeAdded = 0.0;

    const point& min = bb.min();

    vector span = bb.span();

    while (volumeAdded < totalVolume_)
    {
        point p =
            min
          + vector
            (
                span.x()*rndGen.scalar01(),
                span.y()*rndGen.scalar01(),
                span.z()*rndGen.scalar01()
            );

        scalar localSize = cvMesh_.cellSizeControl().cellSize(p);

        scalar localDensity = 1/pow3(max(localSize, VSMALL));

        // Accept possible placements proportional to the relative local density
        if (localDensity/maxDensity_ > rndGen.scalar01())
        {
            // Determine if the point is "wellInside" the domain
            if
            (
                cvMesh_.geometryToConformTo().wellInside
                (
                    p,
                    minimumSurfaceDistanceCoeffSqr_*sqr(localSize)
                )
            )
            {
                initialPoints.push_back(Vb::Point(p.x(), p.y(), p.z()));

                volumeAdded += 1/localDensity;
            }
        }
    }

    return initialPoints;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
