/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

#include "coalParcel.H"

// Kinematic
#include "makeParcelDispersionModels.H"
#include "makeParcelDragModels.H"
#include "makeReactingMultiphaseParcelInjectionModels.H" // MP variant
#include "makeParcelCollisionModels.H"
#include "makeParcelPatchInteractionModels.H"
#include "makeParcelPostProcessingModels.H"

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
    // Kinematic sub-models
    makeParcelDispersionModels(coalParcel);
    makeParcelDragModels(coalParcel);
    makeReactingMultiphaseParcelInjectionModels(coalParcel);
    makeParcelCollisionModels(coalParcel);
    makeParcelPatchInteractionModels(coalParcel);
    makeParcelPostProcessingModels(coalParcel);

    // Thermo sub-models
    makeParcelHeatTransferModels(coalParcel);

    // Reacting sub-models
    makeReactingMultiphaseParcelCompositionModels(coalParcel);
    makeReactingParcelPhaseChangeModels(coalParcel);

    // Reacting multiphase sub-models
    makeReactingMultiphaseParcelDevolatilisationModels(coalParcel);
    makeReactingParcelSurfaceFilmModels(coalParcel);
    makeCoalParcelSurfaceReactionModels(coalParcel);
}


// ************************************************************************* //
