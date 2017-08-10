/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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

#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoPhysicsTypes.H"

#include "rhoConst.H"
#include "perfectFluid.H"
#include "Boussinesq.H"

#include "pureMixture.H"
#include "multiComponentMixture.H"
#include "reactingMixture.H"
#include "SpecieMixture.H"

#include "rhoThermo.H"
#include "rhoReactionThermo.H"
#include "heRhoThermo.H"

#include "solidThermo.H"
#include "heSolidThermo.H"
#include "solidThermoPhysicsTypes.H"

#include "kineticGasEvaporation.H"
#include "constantMelting.H"
#include "constantSolidification.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef
        constTransport
        <
            species::thermo
            <
                hConstThermo
                <
                    rhoConst<specie>
                >,
                sensibleEnthalpy
            >
        > constFluidHThermoPhysics;


    typedef
        constTransport
        <
            species::thermo
            <
                hConstThermo
                <
                    Boussinesq<specie>
                >,
                sensibleEnthalpy
            >
        > BoussinesqFluidEThermoPhysics;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    using namespace meltingEvaporationModels;

    //NOTE: First thermo (dispersed) and second otherThermo (continous)
    // in the phaseProperties: (dispersed to continous)

    // kineticGasEvaporation model definitions
/*
        // multi-component dispersed phase and a pure continous phase
        makeInterfaceDispSpecieMixtureType
        (
            kineticGasEvaporation,
            heRhoThermo,
            rhoReactionThermo,
            multiComponentMixture,
            constIncompressibleGasHThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constFluidHThermoPhysics
        );
*/

        // pure dispersed phase to a multi-component continous phase
        makeInterfaceContSpecieMixtureType
        (
            kineticGasEvaporation,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constFluidHThermoPhysics,
            heRhoThermo,
            rhoReactionThermo,
            multiComponentMixture,
            constIncompressibleGasHThermoPhysics
        );


        // pure dispersed phase and pure continous phase with incompressible gas
        makeInterfacePureType
        (
            kineticGasEvaporation,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constIncompressibleGasHThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constFluidHThermoPhysics
        );

        // pure dispersed phase and pure continous phase with rhoConst gas
        makeInterfacePureType
        (
            kineticGasEvaporation,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constFluidHThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constFluidHThermoPhysics
        );


    // constantMelting model definitions

        // pure dispersed phase and a pure continous phase
        makeInterfacePureType
        (
            constantMelting,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constFluidHThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constFluidHThermoPhysics
        );

        makeInterfacePureType
        (
            constantMelting,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constFluidHThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            BoussinesqFluidEThermoPhysics
        );


        makeInterfacePureType
        (
            constantMelting,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hConstSolidThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            BoussinesqFluidEThermoPhysics
        );

        makeInterfacePureType
        (
            constantMelting,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constFluidHThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constIncompressibleGasHThermoPhysics
        );

        makeInterfaceContSpecieMixtureType
        (
            constantMelting,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constFluidHThermoPhysics,
            heRhoThermo,
            rhoReactionThermo,
            multiComponentMixture,
            constIncompressibleGasHThermoPhysics
        );

/*
        makeInterfaceDispSpecieMixtureType
        (
            constantMelting,
            heRhoThermo,
            rhoReactionThermo,
            multiComponentMixture,
            constIncompressibleGasHThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constFluidHThermoPhysics
        );
*/
    // constantSolidification model definitions

        // pure dispersed phase and a pure continous phase
        makeInterfacePureType
        (
            constantSolidification,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constFluidHThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constFluidHThermoPhysics
        );


        makeInterfacePureType
        (
            constantSolidification,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            BoussinesqFluidEThermoPhysics,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hConstSolidThermoPhysics
        );

/*
    // saturatedEvaporation model definitions

        // multi-component dispersed phase and a pure continous phase
        makeInterfaceDispSpecieMixtureType
        (
            saturatedEvaporation,
            heRhoThermo,
            rhoReactionThermo,
            multiComponentMixture,
            constIncompressibleGasHThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constFluidHThermoPhysics
        );

        // pure dispersed phase and a multi-component continous phase
        makeInterfaceContSpecieMixtureType
        (
            saturatedEvaporation,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constFluidHThermoPhysics,
            heRhoThermo,
            rhoReactionThermo,
            multiComponentMixture,
            constIncompressibleGasHThermoPhysics
        );

        // multi-component dispersed phase and a multi-componen continous phase
        makeSpecieInterfaceSpecieMixtures
        (
            saturatedEvaporation,
            heRhoThermo,
            rhoReactionThermo,
            multiComponentMixture,
            constIncompressibleGasHThermoPhysics,
            heRhoThermo,
            rhoReactionThermo,
            multiComponentMixture,
            constIncompressibleGasHThermoPhysics
        );

        // pure dispersed phase and pure continous phase
        makeInterfacePureType
        (
            saturatedEvaporation,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constIncompressibleGasHThermoPhysics,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constIncompressibleGasHThermoPhysics
        );
*/
}

// ************************************************************************* //
