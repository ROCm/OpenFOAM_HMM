/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "layerParameters.H"
#include "polyBoundaryMesh.H"
#include "unitConversion.H"
#include "refinementSurfaces.H"
#include "searchableSurfaces.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::layerParameters::defaultConcaveAngle = 90;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::layerParameters::layerParameters
(
    const dictionary& dict,
    const polyBoundaryMesh& boundaryMesh
)
:
    numLayers_(boundaryMesh.size(), -1),
    expansionRatio_
    (
        boundaryMesh.size(),
        readScalar(dict.lookup("expansionRatio"))
    ),
    relativeSizes_(dict.lookup("relativeSizes")),
    finalLayerThickness_
    (
        boundaryMesh.size(),
        readScalar(dict.lookup("finalLayerThickness"))
    ),
    minThickness_
    (
        boundaryMesh.size(),
        readScalar(dict.lookup("minThickness"))
    ),
    featureAngle_(readScalar(dict.lookup("featureAngle"))),
    slipFeatureAngle_
    (
        dict.found("slipFeatureAngle")
      ? readScalar(dict.lookup("slipFeatureAngle"))
      : 0.5*featureAngle_
    ),
    concaveAngle_
    (
        dict.lookupOrDefault("concaveAngle", defaultConcaveAngle)
    ),
    nGrow_(readLabel(dict.lookup("nGrow"))),
    nSmoothSurfaceNormals_
    (
        readLabel(dict.lookup("nSmoothSurfaceNormals"))
    ),
    nSmoothNormals_(readLabel(dict.lookup("nSmoothNormals"))),
    nSmoothThickness_(readLabel(dict.lookup("nSmoothThickness"))),
    maxFaceThicknessRatio_
    (
        readScalar(dict.lookup("maxFaceThicknessRatio"))
    ),
    layerTerminationCos_
    (
        Foam::cos(degToRad(0.5*featureAngle_))
    ),
    maxThicknessToMedialRatio_
    (
        readScalar(dict.lookup("maxThicknessToMedialRatio"))
    ),
    minMedianAxisAngleCos_
    (
        Foam::cos(degToRad(readScalar(dict.lookup("minMedianAxisAngle"))))
    ),
    nBufferCellsNoExtrude_
    (
        readLabel(dict.lookup("nBufferCellsNoExtrude"))
    ),
    nSnap_(readLabel(dict.lookup("nRelaxIter"))),
    nLayerIter_(readLabel(dict.lookup("nLayerIter"))),
    nRelaxedIter_(labelMax),
    additionalReporting_(dict.lookupOrDefault("additionalReporting", false))
{
    if (nGrow_ > 0)
    {
        WarningIn("layerParameters::layerParameters(..)")
            << "The nGrow parameter effect has changed with respect to 1.6.x."
            << endl
            << "Please set nGrow=0 for 1.6.x behaviour."
            << endl;
    }

    dict.readIfPresent("nRelaxedIter", nRelaxedIter_);

    if (nLayerIter_ < 0 || nRelaxedIter_ < 0)
    {
        FatalErrorIn("layerParameters::layerParameters(..)")
            << "Layer iterations should be >= 0." << endl
            << "nLayerIter:" << nLayerIter_
            << " nRelaxedIter:" << nRelaxedIter_
            << exit(FatalError);
    }


    const dictionary& layersDict = dict.subDict("layers");

    forAllConstIter(dictionary, layersDict, iter)
    {
        if (iter().isDict())
        {
            const word& key = iter().keyword();
            const labelHashSet patchIDs
            (
                boundaryMesh.patchSet(List<wordRe>(1, key))
            );

            if (patchIDs.size() == 0)
            {
                IOWarningIn("layerParameters::layerParameters(..)", layersDict)
                    << "Layer specification for " << key
                    << " does not match any patch." << endl
                    << "Valid patches are " << boundaryMesh.names() << endl;
            }
            else
            {
                const dictionary& layerDict = iter().dict();

                forAllConstIter(labelHashSet, patchIDs, patchIter)
                {
                    label patchI = patchIter.key();

                    numLayers_[patchI] =
                        readLabel(layerDict.lookup("nSurfaceLayers"));

                    layerDict.readIfPresent
                    (
                        "expansionRatio",
                        expansionRatio_[patchI]
                    );
                    layerDict.readIfPresent
                    (
                        "finalLayerThickness",
                        finalLayerThickness_[patchI]
                    );
                    layerDict.readIfPresent
                    (
                        "minThickness",
                        minThickness_[patchI]
                    );
                }
            }
        }
    }
}


// ************************************************************************* //
