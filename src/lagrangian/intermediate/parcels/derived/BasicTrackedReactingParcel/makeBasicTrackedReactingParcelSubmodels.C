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

#include "BasicTrackedReactingParcel.H"

// Kinematic
#include "makeReactingParcelDispersionModels.H"
#include "makeReactingParcelDragModels.H"
#include "makeReactingParcelInjectionModels.H"
#include "makeReactingParcelPatchInteractionModels.H"
#include "makeReactingParcelPostProcessingModels.H"

// Thermo
#include "makeReactingParcelHeatTransferModels.H"

// Reacting
#include "makeReactingParcelCompositionModels.H"
#include "makeReactingParcelPhaseChangeModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Kinematic sub-models
    makeReactingDispersionModels(BasicTrackedReactingParcel);
    makeReactingDragModels(BasicTrackedReactingParcel);
    makeReactingInjectionModels(BasicTrackedReactingParcel);
    makeReactingPatchInteractionModels(BasicTrackedReactingParcel);
    makeReactingPostProcessingModels(BasicTrackedReactingParcel);

    // Thermo sub-models
    makeReactingHeatTransferModels(BasicTrackedReactingParcel);

    // Reacting sub-models
    makeReactingCompositionModels(BasicTrackedReactingParcel);
    makeReactingPhaseChangeModels(BasicTrackedReactingParcel);
};


// ************************************************************************* //
