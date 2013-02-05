/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "automatic.H"
#include "addToRunTimeSelectionTable.H"
#include "triSurfaceMesh.H"
#include "vtkSurfaceWriter.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(automatic, 0);
    addToRunTimeSelectionTable(cellSizeCalculationType, automatic, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::automatic::smoothField(triSurfaceScalarField& surf)
{
    label nSmoothingIterations = 10;

    for (label iter = 0; iter < nSmoothingIterations; ++iter)
    {
        const pointField& faceCentres = surface_.faceCentres();

        forAll(surf, sI)
        {
            const labelList& faceFaces = surface_.faceFaces()[sI];

            const point& fC = faceCentres[sI];
            const scalar value = surf[sI];

            scalar newValue = 0;
            scalar totalDist = 0;

            label nFaces = 0;

            forAll(faceFaces, fI)
            {
                const label faceLabel = faceFaces[fI];
                const point& faceCentre = faceCentres[faceLabel];

                const scalar faceValue = surf[faceLabel];
                const scalar distance = mag(faceCentre - fC);

                newValue += faceValue/distance;

                totalDist += 1.0/distance;

                if (value < faceValue)
                {
                    nFaces++;
                }
            }

            // Do not smooth out the peak values
            if (nFaces == faceFaces.size())
            {
                //continue;
            }

            surf[sI] = newValue/totalDist;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::automatic::automatic
(
    const dictionary& cellSizeCalcTypeDict,
    const triSurfaceMesh& surface,
    const scalar& defaultCellSize
)
:
    cellSizeCalculationType
    (
        typeName,
        cellSizeCalcTypeDict,
        surface,
        defaultCellSize
    ),
    coeffsDict_(cellSizeCalcTypeDict.subDict(typeName + "Coeffs")),
    surface_(surface),
    surfaceName_(surface.searchableSurface::name()),
    readCurvature_(Switch(coeffsDict_.lookup("curvature"))),
    curvatureFile_(coeffsDict_.lookup("curvatureFile")),
    readFeatureProximity_(Switch(coeffsDict_.lookup("featureProximity"))),
    featureProximityFile_(coeffsDict_.lookup("featureProximityFile")),
    readInternalCloseness_(Switch(coeffsDict_.lookup("internalCloseness"))),
    internalClosenessFile_(coeffsDict_.lookup("internalClosenessFile")),
    curvatureCellSizeCoeff_
    (
        readScalar(coeffsDict_.lookup("curvatureCellSizeCoeff"))
    ),
    maximumCellSize_
    (
        readScalar(coeffsDict_.lookup("maximumCellSizeCoeff"))*defaultCellSize_
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::triSurfaceScalarField Foam::automatic::load()
{
    Info<< indent << "Calculating cell size on surface: "
        << surfaceName_ << endl;

    triSurfaceScalarField surfaceCellSize
    (
        IOobject
        (
            surfaceName_ + ".cellSize",
            surface_.searchableSurface::time().constant(),
            "triSurface",
            surface_.searchableSurface::time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        surface_,
        dimLength,
        scalarField(surface_.size(), maximumCellSize_)
    );

    if (readCurvature_)
    {
        Info<< indent << "Reading curvature         : "
            << curvatureFile_ << endl;

        triSurfacePointScalarField curvature
        (
            IOobject
            (
                curvatureFile_,
                surface_.searchableSurface::time().constant(),
                "triSurface",
                surface_.searchableSurface::time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            surface_,
            dimLength,
            true
        );

        const List<labelledTri>& localFaces = surface_.localFaces();
        const labelList& meshPoints = surface_.meshPoints();

        forAll(surfaceCellSize, fI)
        {
            const labelList& facePoints = localFaces[fI].triFaceFace();

            scalar interpolatedCurvatureToFace = 0.0;

            forAll(facePoints, fpI)
            {
                interpolatedCurvatureToFace
                    += curvature[meshPoints[facePoints[fpI]]];
            }

            interpolatedCurvatureToFace /= facePoints.size();

            surfaceCellSize[fI] =
                min
                (
                    1.0
                   /max
                    (
                        (1.0/curvatureCellSizeCoeff_)
                       *interpolatedCurvatureToFace,
                        1.0/maximumCellSize_
                    ),
                    surfaceCellSize[fI]
                );
        }
    }

    if (readInternalCloseness_)
    {
        Info<< indent << "Reading internal closeness: "
            << internalClosenessFile_ << endl;

        triSurfaceScalarField internalCloseness
        (
            IOobject
            (
                internalClosenessFile_,
                surface_.searchableSurface::time().constant(),
                "triSurface",
                surface_.searchableSurface::time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            surface_,
            dimLength,
            true
        );

        forAll(surfaceCellSize, fI)
        {
            surfaceCellSize[fI] =
                min
                (
                    internalCloseness[fI],
                    surfaceCellSize[fI]
                );
        }
    }

    if (readFeatureProximity_)
    {
        Info<< indent << "Reading feature proximity : "
            << featureProximityFile_ << endl;

        triSurfaceScalarField featureProximity
        (
            IOobject
            (
                featureProximityFile_,
                surface_.searchableSurface::time().constant(),
                "triSurface",
                surface_.searchableSurface::time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            surface_,
            dimLength,
            true
        );

        forAll(surfaceCellSize, fI)
        {
            surfaceCellSize[fI] =
                min
                (
                    featureProximity[fI],
                    surfaceCellSize[fI]
                );
        }
    }

    smoothField(surfaceCellSize);

    surfaceCellSize.write();

    if (debug)
    {
        faceList faces(surface_.size());

        forAll(surface_, fI)
        {
            faces[fI] = surface_.triSurface::operator[](fI).triFaceFace();
        }

        vtkSurfaceWriter().write
        (
            surface_.searchableSurface::time().constant()/"triSurface",
            surfaceName_,
            surface_.points(),
            faces,
            "cellSize",
            surfaceCellSize,
            false,
            true
        );
    }

    return surfaceCellSize;
}


// ************************************************************************* //
