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

#include "makeReactionThermo.H"

#include "psiReactionThermo.H"
#include "hePsiReactionThermo.H"

#include "perfectGas.H"

#include "hConstThermo.H"
#include "janafThermo.H"
#include "sensibleEnthalpy.H"
#include "specieThermo.H"

#include "constTransport.H"
#include "sutherlandTransport.H"

#include "homogeneousMixture.H"
#include "inhomogeneousMixture.H"
#include "veryInhomogeneousMixture.H"
#include "dieselMixture.H"
#include "multiComponentMixture.H"
#include "reactingMixture.H"
#include "singleStepReactingMixture.H"

#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// constTransport, hConstThermo

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiReactionThermo,
    homogeneousMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiReactionThermo,
    inhomogeneousMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiReactionThermo,
    veryInhomogeneousMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas
);


// sutherlandTransport, hConstThermo

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiReactionThermo,
    homogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiReactionThermo,
    inhomogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiReactionThermo,
    veryInhomogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas
);


// sutherlandTransport, janafThermo

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiReactionThermo,
    homogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiReactionThermo,
    inhomogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiReactionThermo,
    veryInhomogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas
);


makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiReactionThermo,
    dieselMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas
);

// Multi-component thermo

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiReactionThermo,
    multiComponentMixture,
    constGasThermoPhysics
);

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiReactionThermo,
    multiComponentMixture,
    gasThermoPhysics
);


// Multi-component reaction thermo

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiReactionThermo,
    reactingMixture,
    constGasThermoPhysics
);

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiReactionThermo,
    reactingMixture,
    gasThermoPhysics
);

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiReactionThermo,
    singleStepReactingMixture,
    gasThermoPhysics
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
