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

// Kinematic
#include "makeThermoParcelForces.H" // thermo variant
#include "makeParcelDispersionModels.H"
#include "makeReactingMultiphaseParcelInjectionModels.H" // MP variant
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
    typedef coalCloud::cloudType coalCloud_R;
    typedef coalCloud_R::cloudType coalCloud_T;
    typedef coalCloud_T::cloudType coalCloud_K;

    // Kinematic sub-models
    makeThermoParcelForces(coalCloud_K);
    makeParcelDispersionModels(coalCloud_K);
    makeReactingMultiphaseParcelInjectionModels(coalCloud_K);
    makeParcelPatchInteractionModels(coalCloud_K);
    makeParcelPostProcessingModels(coalCloud_K);

    // Thermo sub-models
    makeParcelHeatTransferModels(coalCloud_T);

    // Reacting sub-models
    makeReactingMultiphaseParcelCompositionModels(coalCloud_R);
    makeReactingParcelPhaseChangeModels(coalCloud_R);

    // Reacting multiphase sub-models
    makeReactingMultiphaseParcelDevolatilisationModels(coalCloud);
    makeReactingParcelSurfaceFilmModels(coalCloud_K);
    makeCoalParcelSurfaceReactionModels(coalCloud);
}


// ************************************************************************* //
