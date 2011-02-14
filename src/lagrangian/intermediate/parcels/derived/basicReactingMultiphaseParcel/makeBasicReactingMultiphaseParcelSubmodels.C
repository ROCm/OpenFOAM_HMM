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

#include "basicReactingMultiphaseParcel.H"

// Kinematic
#include "makeParcelDispersionModels.H"
#include "makeReactingMultiphaseParcelInjectionModels.H" // MP variant
#include "makeParcelCollisionModels.H"
#include "makeParcelPatchInteractionModels.H"
#include "makeParcelPostProcessingModels.H"

// Thermodynamic
#include "makeParcelHeatTransferModels.H"

// Reacting
#include "makeReactingMultiphaseParcelCompositionModels.H" // MP Variant
#include "makeReactingParcelPhaseChangeModels.H"
#include "makeReactingParcelSurfaceFilmModels.H"

// Reacting multiphase
#include "makeReactingMultiphaseParcelDevolatilisationModels.H"
#include "makeReactingMultiphaseParcelSurfaceReactionModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Kinematic sub-models
    makeParcelDispersionModels(basicReactingMultiphaseParcel);
    makeReactingMultiphaseParcelInjectionModels(basicReactingMultiphaseParcel);
    makeParcelCollisionModels(basicReactingMultiphaseParcel);
    makeParcelPatchInteractionModels(basicReactingMultiphaseParcel);
    makeParcelPostProcessingModels(basicReactingMultiphaseParcel);

    // Thermo sub-models
    makeParcelHeatTransferModels(basicReactingMultiphaseParcel);

    // Reacting sub-models
    makeReactingMultiphaseParcelCompositionModels
    (
        basicReactingMultiphaseParcel
    );
    makeReactingParcelPhaseChangeModels(basicReactingMultiphaseParcel);

    // Reacting multiphase sub-models
    makeReactingMultiphaseParcelDevolatilisationModels
    (
        basicReactingMultiphaseParcel
    );
    makeReactingParcelSurfaceFilmModels
    (
        basicReactingMultiphaseParcel
    );
    makeReactingMultiphaseParcelSurfaceReactionModels
    (
        basicReactingMultiphaseParcel
    );
}


// ************************************************************************* //
