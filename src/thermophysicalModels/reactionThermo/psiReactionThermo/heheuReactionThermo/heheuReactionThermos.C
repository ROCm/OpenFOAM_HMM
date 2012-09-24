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

#include "psiuReactionThermo.H"
#include "heheuReactionThermo.H"

#include "makeReactionThermo.H"
#include "addToRunTimeSelectionTable.H"

#include "perfectGas.H"

#include "hConstThermo.H"
#include "janafThermo.H"
#include "absoluteEnthalpy.H"
#include "thermo.H"

#include "constTransport.H"
#include "sutherlandTransport.H"

#include "homogeneousMixture.H"
#include "inhomogeneousMixture.H"
#include "veryInhomogeneousMixture.H"
#include "multiComponentMixture.H"
#include "egrMixture.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * h-hu-Thermos * * * * * * * * * * * * * * * //

makeReactionThermo
(
    psiThermo,
    psiuReactionThermo,
    heheuReactionThermo,
    homogeneousMixture,
    constTransport,
    absoluteEnthalpy,
    hConstThermo,
    perfectGas
);

makeReactionThermo
(
    psiThermo,
    psiuReactionThermo,
    heheuReactionThermo,
    inhomogeneousMixture,
    constTransport,
    absoluteEnthalpy,
    hConstThermo,
    perfectGas
);

makeReactionThermo
(
    psiThermo,
    psiuReactionThermo,
    heheuReactionThermo,
    veryInhomogeneousMixture,
    constTransport,
    absoluteEnthalpy,
    hConstThermo,
    perfectGas
);

makeReactionThermo
(
    psiThermo,
    psiuReactionThermo,
    heheuReactionThermo,
    homogeneousMixture,
    sutherlandTransport,
    absoluteEnthalpy,
    janafThermo,
    perfectGas
);

makeReactionThermo
(
    psiThermo,
    psiuReactionThermo,
    heheuReactionThermo,
    inhomogeneousMixture,
    sutherlandTransport,
    absoluteEnthalpy,
    janafThermo,
    perfectGas
);

makeReactionThermo
(
    psiThermo,
    psiuReactionThermo,
    heheuReactionThermo,
    veryInhomogeneousMixture,
    sutherlandTransport,
    absoluteEnthalpy,
    janafThermo,
    perfectGas
);

makeReactionThermo
(
    psiThermo,
    psiuReactionThermo,
    heheuReactionThermo,
    egrMixture,
    constTransport,
    absoluteEnthalpy,
    hConstThermo,
    perfectGas
);

makeReactionThermo
(
    psiThermo,
    psiuReactionThermo,
    heheuReactionThermo,
    egrMixture,
    sutherlandTransport,
    absoluteEnthalpy,
    janafThermo,
    perfectGas
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
