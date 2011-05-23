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

#include "BasicSprayParcel.H"

// Kinematic
#include "makeSprayParcelDispersionModels.H"
#include "makeSprayParcelDragModels.H"
#include "makeSprayParcelInjectionModels.H"
#include "makeSprayParcelPatchInteractionModels.H"
#include "makeSprayParcelPostProcessingModels.H"

// Thermodynamic
#include "makeSprayParcelHeatTransferModels.H"

// Reacting
#include "makeSprayParcelCompositionModels.H"
#include "makeSprayParcelPhaseChangeModels.H"

// Spray
#include "makeSprayParcelAtomizationModels.H"
#include "makeSprayParcelBreakupModels.H"
#include "makeSprayParcelCollisionModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Kinematic sub-models
    makeSprayDispersionModels(BasicSprayParcel);
    makeSprayDragModels(BasicSprayParcel);
    makeSprayInjectionModels(BasicSprayParcel);
    makeSprayPatchInteractionModels(BasicSprayParcel);
    makeSprayPostProcessingModels(BasicSprayParcel);

    // Thermo sub-models
    makeSprayHeatTransferModels(BasicSprayParcel);

    // Reacting sub-models
    makeSprayCompositionModels(BasicSprayParcel);
    makeSprayPhaseChangeModels(BasicSprayParcel);

    // Spray sub-models
    makeSprayAtomizationModels(BasicSprayParcel);
    makeSprayBreakupModels(BasicSprayParcel);
    makeSprayCollisionModels(BasicSprayParcel);

};


// ************************************************************************* //
