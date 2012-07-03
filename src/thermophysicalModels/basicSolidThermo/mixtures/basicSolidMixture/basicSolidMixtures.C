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

Description
    Mixture instantiation

\*---------------------------------------------------------------------------*/


#include "basicMixture.H"
#include "makeBasicMixture.H"


#include "constRho.H"

#include "constSolidThermo.H"
#include "exponentialSolidThermo.H"

#include "constIsoSolidTransport.H"
#include "constAnIsoSolidTransport.H"
#include "exponentialSolidTransport.H"

#include "constSolidRad.H"

#include "sensibleInternalEnergy.H"
#include "sensibleEnthalpy.H"
#include "specieThermo.H"

#include "pureSolidMixture.H"
#include "multiComponentSolidMixture.H"
#include "reactingSolidMixture.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

makeBasicMixture
(
    pureSolidMixture,
    constIsoSolidTransport,
    sensibleEnthalpy,
    constSolidThermo,
    constRho
);

makeBasicMixture
(
    pureSolidMixture,
    constAnIsoSolidTransport,
    sensibleEnthalpy,
    constSolidThermo,
    constRho
);

makeBasicMixture
(
    pureSolidMixture,
    exponentialSolidTransport,
    sensibleEnthalpy,
    exponentialSolidThermo,
    constRho
);

makeBasicMixture
(
    multiComponentSolidMixture,
    constIsoSolidTransport,
    sensibleEnthalpy,
    constSolidThermo,
    constRho
);

makeBasicMixture
(
    reactingSolidMixture,
    constIsoSolidTransport,
    sensibleEnthalpy,
    constSolidThermo,
    constRho
);


/*
makeBasicMixture
(
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    eConstThermo,
    perfectGas
);

makeBasicMixture
(
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hConstThermo,
    perfectGas
);

makeBasicMixture
(
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    hConstThermo,
    perfectGas
);

makeBasicMixture
(
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    janafThermo,
    perfectGas
);

makeBasicMixture
(
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas
);

makeBasicMixture
(
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas
);

makeBasicMixture
(
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas
);

makeBasicMixture
(
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    incompressible
);

makeBasicPolyMixture
(
    pureMixture,
    3,
    sensibleEnthalpy
);

makeBasicPolyMixture
(
    pureMixture,
    8,
    sensibleEnthalpy
);


makeBasicMixture
(
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    isobaricPerfectGas
);

makeBasicMixture
(
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    isobaricPerfectGas
);

makeBasicMixture
(
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    isobaricPerfectGas
);
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
