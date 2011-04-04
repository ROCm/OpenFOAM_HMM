/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
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

#include "pointFile.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pointFile, 0);
addToRunTimeSelectionTable(initialPointsMethod, pointFile, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pointFile::pointFile
(
    const dictionary& initialPointsDict,
    const conformalVoronoiMesh& cvMesh
)
:
    initialPointsMethod(typeName, initialPointsDict, cvMesh),
    pointFileName_(detailsDict().lookup("pointFile"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

std::vector<Vb::Point> pointFile::initialPoints() const
{
    pointIOField points
    (
        IOobject
        (
            pointFileName_.name(),
            pointFileName_.path(),
            cvMesh_.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info<< "    Inserting points from file " << pointFileName_ << endl;

    if (points.empty())
    {
        FatalErrorIn("std::vector<Vb::Point> pointFile::initialPoints() const")
            << "Point file contain no points"
            << exit(FatalError) << endl;
    }

    // Filter the points to be only those on this processor
    boolList procPt(points.size(), false);

    forAll(points, ptI)
    {
        procPt[ptI] = cvMesh_.geometryToConformTo().positionOnThisProc
        (
            points[ptI]
        );
    }

    if (Pstream::parRun())
    {
        List<boolList> allProcPt(Pstream::nProcs());

        allProcPt[Pstream::myProcNo()] = procPt;

        Pstream::gatherList(allProcPt);

        Pstream::scatterList(allProcPt);

        // forAll(procPt, ptI)
        // {
        //     label nProcsTrue = 0;

        //     forAll(allProcPt, procI)
        //     {
        //         if (allProcPt[procI][ptI])
        //         {
        //             nProcsTrue++;
        //         }
        //     }

        //     if (nProcsTrue != 1)
        //     {
        //         Info<< points[ptI] << " " << nProcsTrue << " true" << endl;

        //         if (nProcsTrue > 1)
        //         {
        //             forAll(allProcPt, procI)
        //             {
        //                 allProcPt[procI][ptI] = false;
        //             }
        //         }
        //     }
        // }

        forAll(procPt, ptI)
        {
            bool foundAlready = false;

            forAll(allProcPt, procI)
            {
                // If a processor with a lower index has found this point to
                // insert already, defer to it and don't insert.
                if (foundAlready)
                {
                    allProcPt[procI][ptI] = false;
                }
                else if (allProcPt[procI][ptI])
                {
                    foundAlready = true;
                }
            }
        }

        procPt = allProcPt[Pstream::myProcNo()];
    }

    inplaceSubset(procPt, points);

    std::vector<Vb::Point> initialPoints;

    Field<bool> insidePoints = cvMesh_.geometryToConformTo().wellInside
    (
        points,
        minimumSurfaceDistanceCoeffSqr_
       *sqr
        (
            cvMesh_.cellSizeControl().cellSize
            (
                points,
                List<bool>(points.size(), false)
            )
        )
    );

    forAll(insidePoints, i)
    {
        if (insidePoints[i])
        {
            const point& p(points[i]);

            initialPoints.push_back(Vb::Point(p.x(), p.y(), p.z()));
        }
    }

    label nPointsRejected = points.size() - initialPoints.size();

    if (Pstream::parRun())
    {
        reduce(nPointsRejected, sumOp<label>());
    }

    if (nPointsRejected)
    {
        Info<< "    " << nPointsRejected << " points rejected from "
            << pointFileName_.name() << endl;
    }

    return initialPoints;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
