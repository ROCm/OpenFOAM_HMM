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

#include "uniformGrid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(uniformGrid, 0);
addToRunTimeSelectionTable(initialPointsMethod, uniformGrid, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

uniformGrid::uniformGrid
(
    const dictionary& initialPointsDict,
    const conformalVoronoiMesh& cvMesh
)
:
    initialPointsMethod(typeName, initialPointsDict, cvMesh),
    initialCellSize_(readScalar(detailsDict().lookup("initialCellSize"))),
    randomiseInitialGrid_(detailsDict().lookup("randomiseInitialGrid")),
    randomPerturbationCoeff_
    (
        readScalar(detailsDict().lookup("randomPerturbationCoeff"))
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

std::vector<Vb::Point> uniformGrid::initialPoints() const
{
    const boundBox& bb = cvMesh_.geometryToConformTo().bounds();

    scalar x0 = bb.min().x();
    scalar xR = bb.max().x() - x0;
    int ni = int(xR/initialCellSize_) + 1;

    scalar y0 = bb.min().y();
    scalar yR = bb.max().y() - y0;
    int nj = int(yR/initialCellSize_) + 1;

    scalar z0 = bb.min().z();
    scalar zR = bb.max().z() - z0;
    int nk = int(zR/initialCellSize_) + 1;

    Info<< "    Is this actually uniform?  or is it fitting the span with an "
        << "integer number?" << endl;

    vector delta(xR/ni, yR/nj, zR/nk);

    delta *= pow((1.0),-(1.0/3.0));

    Random rndGen(1735621);
    scalar pert = randomPerturbationCoeff_*cmptMin(delta);

    std::vector<Vb::Point> initialPoints;

    for (int i = 0; i < ni; i++)
    {
        for (int j = 0; j < nj; j++)
        {
            // Generating, testing and adding points one line at a time to
            // reduce the memory requirement for cases with bounding boxes that
            // are very large in comparison to the volume to be filled

            label pI = 0;

            pointField points(nk);

            for (int k = 0; k < nk; k++)
            {
                point p
                (
                    x0 + i*delta.x(),
                    y0 + j*delta.y(),
                    z0 + k*delta.z()
                );

                if (randomiseInitialGrid_)
                {
                    p.x() += pert*(rndGen.scalar01() - 0.5);
                    p.y() += pert*(rndGen.scalar01() - 0.5);
                    p.z() += pert*(rndGen.scalar01() - 0.5);
                }

                points[pI++] = p;
            }

            Field<bool> insidePoints = cvMesh_.geometryToConformTo().wellInside
            (
                points,
                minimumSurfaceDistance_*minimumSurfaceDistance_
            );

            forAll(insidePoints, i)
            {
                if (insidePoints[i])
                {
                    const point& p(points[i]);

                    initialPoints.push_back(Vb::Point(p.x(), p.y(), p.z()));
                }
            }
        }
    }

    return initialPoints;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
