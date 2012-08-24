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


#include "incompressible.H"

#include "hConstThermo.H"
#include "hExponentialThermo.H"

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
    hConstThermo,
    incompressible
);

makeBasicMixture
(
    pureSolidMixture,
    constAnIsoSolidTransport,
    sensibleEnthalpy,
    hConstThermo,
    incompressible
);

makeBasicMixture
(
    pureSolidMixture,
    exponentialSolidTransport,
    sensibleEnthalpy,
    hExponentialThermo,
    incompressible
);

makeBasicMixture
(
    multiComponentSolidMixture,
    constIsoSolidTransport,
    sensibleEnthalpy,
    hConstThermo,
    incompressible
);

makeBasicMixture
(
    reactingSolidMixture,
    constIsoSolidTransport,
    sensibleEnthalpy,
    hConstThermo,
    incompressible
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
