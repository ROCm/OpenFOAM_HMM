/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011 OpenCFD Ltd.
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

#include "coalCloud.H"

#include "makeParcelCloudFunctionObjects.H"

// Kinematic
#include "makeThermoParcelForces.H" // thermo variant
#include "makeParcelDispersionModels.H"
#include "makeReactingMultiphaseParcelInjectionModels.H" // MP variant
#include "makeParcelPatchInteractionModels.H"

// Thermodynamic
#include "makeParcelHeatTransferModels.H"

// Reacting
#include "makeReactingMultiphaseParcelCompositionModels.H" // MP variant
#include "makeReactingParcelPhaseChangeModels.H"
#include "makeReactingParcelSurfaceFilmModels.H"

// Reacting multiphase
#include "makeReactingMultiphaseParcelDevolatilisationModels.H"

// Coal specific
#include "makeCoalParcelSurfaceReactionModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeParcelCloudFunctionObjects(coalCloud);

    // Kinematic sub-models
    makeThermoParcelForces(coalCloud);
    makeParcelDispersionModels(coalCloud);
    makeReactingMultiphaseParcelInjectionModels(coalCloud);
    makeParcelPatchInteractionModels(coalCloud);

    // Thermo sub-models
    makeParcelHeatTransferModels(coalCloud);

    // Reacting sub-models
    makeReactingMultiphaseParcelCompositionModels(coalCloud);
    makeReactingParcelPhaseChangeModels(coalCloud);

    // Reacting multiphase sub-models
    makeReactingMultiphaseParcelDevolatilisationModels(coalCloud);
    makeReactingParcelSurfaceFilmModels(coalCloud);
    makeCoalParcelSurfaceReactionModels(coalCloud);
}


// ************************************************************************* //
