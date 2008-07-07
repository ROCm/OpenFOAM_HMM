/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

#include "layerParameters.H"
#include "polyBoundaryMesh.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::layerParameters::defaultConcaveAngle = 90;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::layerParameters::layerParameters
(
    const dictionary& dict,
    const labelList& numLayers
)
:
    numLayers_(numLayers),
    expansionRatio_
    (
        numLayers_.size(),
        readScalar(dict.lookup("expansionRatio"))
    ),
    finalLayerRatio_
    (
        numLayers_.size(),
        readScalar(dict.lookup("finalLayerRatio"))
    ),
    minThickness_
    (
        numLayers_.size(),
        readScalar(dict.lookup("minThickness"))
    ),
    featureAngle_(readScalar(dict.lookup("featureAngle"))),
    concaveAngle_
    (
        dict.found("concaveAngle")
      ? readScalar(dict.lookup("concaveAngle"))
      : defaultConcaveAngle
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
        Foam::cos
        (
            0.5
          * featureAngle_
          * mathematicalConstant::pi/180.
        )
    ),
    maxThicknessToMedialRatio_
    (
        readScalar(dict.lookup("maxThicknessToMedialRatio"))
    ),
    minMedianAxisAngleCos_
    (
        Foam::cos(readScalar(dict.lookup("minMedianAxisAngle")))
      * mathematicalConstant::pi/180.
    ),
    nBufferCellsNoExtrude_
    (
        readLabel(dict.lookup("nBufferCellsNoExtrude"))
    ),
    nSnap_(readLabel(dict.lookup("nSnap")))
{}


// Construct from dictionary
Foam::layerParameters::layerParameters
(
    const dictionary& dict,
    const polyBoundaryMesh& boundaryMesh
)
:
    numLayers_(boundaryMesh.size(), 0),
    expansionRatio_
    (
        boundaryMesh.size(),
        readScalar(dict.lookup("expansionRatio"))
    ),
    finalLayerRatio_
    (
        boundaryMesh.size(),
        readScalar(dict.lookup("finalLayerRatio"))
    ),
    minThickness_
    (
        boundaryMesh.size(),
        readScalar(dict.lookup("minThickness"))
    ),
    featureAngle_(readScalar(dict.lookup("featureAngle"))),
    concaveAngle_
    (
        dict.found("concaveAngle")
      ? readScalar(dict.lookup("concaveAngle"))
      : defaultConcaveAngle
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
        Foam::cos
        (
            0.5
          * featureAngle_
          * mathematicalConstant::pi/180.
        )
    ),
    maxThicknessToMedialRatio_
    (
        readScalar(dict.lookup("maxThicknessToMedialRatio"))
    ),
    minMedianAxisAngleCos_
    (
        Foam::cos(readScalar(dict.lookup("minMedianAxisAngle")))
      * mathematicalConstant::pi/180.
    ),
    nBufferCellsNoExtrude_
    (
        readLabel(dict.lookup("nBufferCellsNoExtrude"))
    ),
    nSnap_(readLabel(dict.lookup("nRelaxIter")))
{
    const dictionary& layersDict = dict.subDict("layers");

    forAllConstIter(dictionary, layersDict, iter)
    {
        const word& key = iter().keyword();

        if (layersDict.isDict(key))
        {
            label patchI = boundaryMesh.findPatchID(key);

            if (patchI == -1)
            {
                FatalErrorIn
                (
                    "layerParameters::layerParameters"
                    "(const dictionary&, const polyBoundaryMesh&)"
                )   << "Specified illegal patch " << key
                    << " in layer dictionary." << endl
                    << "Valid patch names are " << boundaryMesh.names()
                    << exit(FatalError);
            }

            const dictionary& layerDict = layersDict.subDict(key);

            numLayers_[patchI] =
                readLabel(layerDict.lookup("nSurfaceLayers"));

            //- Patch-wise layer parameters disabled for now. Just remove
            //  settings in initialiser list and uncomment below.
            //expansionRatio_[patchI] =
            //    readScalar(layerDict.lookup("expansionRatio"));
            //finalLayerRatio_[patchI] =
            //    readScalar(layerDict.lookup("finalLayerRatio"));
            //minThickness_[patchI] =
            //    readScalar(layerDict.lookup("minThickness"));
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
