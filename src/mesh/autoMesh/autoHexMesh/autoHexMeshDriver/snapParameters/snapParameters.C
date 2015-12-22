/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "snapParameters.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::snapParameters::snapParameters(const dictionary& dict)
:
    nSmoothPatch_(readLabel(dict.lookup("nSmoothPatch"))),
    nSmoothInternal_(dict.lookupOrDefault("nSmoothInternal", 0)),
    snapTol_(readScalar(dict.lookup("tolerance"))),
    nSmoothDispl_(readLabel(dict.lookup("nSolveIter"))),
    nSnap_(readLabel(dict.lookup("nRelaxIter"))),
    nFeatureSnap_(dict.lookupOrDefault("nFeatureSnapIter", -1)),
    explicitFeatureSnap_(dict.lookupOrDefault("explicitFeatureSnap", true)),
    implicitFeatureSnap_(dict.lookupOrDefault("implicitFeatureSnap", false)),
    multiRegionFeatureSnap_
    (
        dict.lookupOrDefault("multiRegionFeatureSnap", false)
    ),
    detectNearSurfacesSnap_
    (
        dict.lookupOrDefault("detectNearSurfacesSnap", true)
    ),
    strictRegionSnap_
    (
        dict.lookupOrDefault("strictRegionSnap", false)
    ),
    detectBaffles_(dict.lookupOrDefault("detectBaffles", true)),
    baffleFeaturePoints_(dict.lookupOrDefault("baffleFeaturePoints", false)),
    releasePoints_(dict.lookupOrDefault("releasePoints", false)),
    stringFeatures_(dict.lookupOrDefault("stringFeatures", true)),
    avoidDiagonal_(dict.lookupOrDefault("avoidDiagonal", false)),
    nFaceSplitInterval_
    (
        dict.lookupOrDefault("nFaceSplitInterval", labelMin)
    ),
    concaveAngle_(dict.lookupOrDefault("concaveAngle", 45)),
    minAreaRatio_(dict.lookupOrDefault("minAreaRatio", 0.3))
{}


// ************************************************************************* //
