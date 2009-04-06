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

    std::vector<Vb::Point> initialPoints;

    forAll(points, i)
    {
        const point& p = points[i];

        // TODO Check if inside the surface

        initialPoints.push_back(Vb::Point(p.x(), p.y(), p.z()));
    }

    label nPointsRejected = points.size() - initialPoints.size();

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
