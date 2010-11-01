/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

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
    minCellSize_
    (
        detailsDict().lookupOrDefault<scalar>("minCellSize", GREAT)
    ),
    minCellSizeLimit_
    (
        detailsDict().lookupOrDefault<scalar>("minCellSizeLimit", 0.0)
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

    label trialPoints = 0;

    scalar maxDensity = 1/pow3(max(minCellSize_, SMALL));

    while (volumeAdded < totalVolume_)
    {
        trialPoints++;

        point p =
            min
          + vector
            (
                span.x()*rndGen.scalar01(),
                span.y()*rndGen.scalar01(),
                span.z()*rndGen.scalar01()
            );

        scalar localSize = cvMesh_.cellSizeControl().cellSize(p);

        if (localSize < minCellSize_)
        {
            minCellSize_ = max(localSize, minCellSizeLimit_);

            // 1/(minimum cell size)^3, gives the maximum permissible point
            // density
            maxDensity = 1/pow3(max(minCellSize_, SMALL));
        }

        scalar localDensity = 1/pow3(max(localSize, SMALL));

        // Accept possible placements proportional to the relative local density
        if (localDensity/maxDensity > rndGen.scalar01())
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

    Info<< nl << "    " << typeName << " - "
        << trialPoints << " locations queried ("
        << scalar(initialPoints.size())/scalar(trialPoints)
        << " success rate).  minCellSize " << minCellSize_
        << endl;

    return initialPoints;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
