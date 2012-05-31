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

InClass
    Foam::basicMixture

Description

\*---------------------------------------------------------------------------*/

#ifndef makeBasicMixture_H
#define makeBasicMixture_H

#include "basicMixture.H"
#include "sensibleEnthalpy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeBasicMixture(Mixture,Transport,Type,Thermo,EqnOfState)            \
                                                                              \
typedef                                                                       \
    Mixture<Transport<specieThermo<Thermo<EqnOfState>, Type> > >              \
    Mixture##Transport##Type##Thermo##EqnOfState;                             \
                                                                              \
defineTemplateTypeNameAndDebugWithName                                        \
    (Mixture##Transport##Type##Thermo##EqnOfState,                            \
    #Mixture"<"#Transport"<specieThermo<"#Thermo"<"#EqnOfState">,"#Type">>>", 0)


#define makeBasicPolyMixture(Mixture,Order,Type)                              \
                                                                              \
typedef polynomialTransport                                                   \
<                                                                             \
    specieThermo                                                              \
    <                                                                         \
        hPolynomialThermo                                                     \
        <                                                                     \
            icoPolynomial<Order>,                                             \
            Order                                                             \
        >,                                                                    \
        Type                                                                  \
    >,                                                                        \
    Order                                                                     \
> icoPoly##Order##Type##ThermoPhysics;                                        \
                                                                              \
typedef Mixture<icoPoly##Order##Type##ThermoPhysics>                          \
    Mixture##icoPoly##Order##Type##ThermoPhysics;                             \
                                                                              \
defineTemplateTypeNameAndDebugWithName                                        \
    (Mixture##icoPoly##Order##Type##ThermoPhysics,                            \
    #Mixture"<icoPoly"#Order#Type"ThermoPhysics>", 0)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
