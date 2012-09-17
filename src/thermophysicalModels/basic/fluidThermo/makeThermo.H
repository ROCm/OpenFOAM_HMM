/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::fluidThermo

Description
    Macros for creating 'basic' density-based thermo packages

\*---------------------------------------------------------------------------*/

#ifndef makeThermo_H
#define makeThermo_H

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeThermo(BaseThermo,Cthermo,Mixture,Transport,Type,Thermo,EqnOfState)\
                                                                              \
typedef Cthermo                                                               \
<                                                                             \
    Mixture                                                                   \
    <                                                                         \
        Transport                                                             \
        <                                                                     \
            specieThermo                                                      \
            <                                                                 \
                Thermo                                                        \
                <                                                             \
                    EqnOfState                                                \
                >,                                                            \
                Type                                                          \
            >                                                                 \
        >                                                                     \
    >                                                                         \
>   Cthermo##Mixture##Transport##Type##Thermo##EqnOfState;                    \
                                                                              \
defineTemplateTypeNameAndDebugWithName                                        \
(                                                                             \
    Cthermo##Mixture##Transport##Type##Thermo##EqnOfState,                    \
    #Cthermo                                                                  \
    "<"                                                                       \
        #Mixture                                                              \
        "<"                                                                   \
            #Transport                                                        \
            "<specieThermo<"                                                  \
                #Thermo                                                       \
                "<"                                                           \
                    #EqnOfState                                               \
                ">,"                                                          \
                #Type                                                         \
        ">>>>",                                                               \
    0                                                                         \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    BaseThermo,                                                               \
    Cthermo##Mixture##Transport##Type##Thermo##EqnOfState,                    \
    fvMesh                                                                    \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    fluidThermo,                                                              \
    Cthermo##Mixture##Transport##Type##Thermo##EqnOfState,                    \
    fvMesh                                                                    \
);


#define makePolyThermo(BaseThermo,Cthermo,Mixture,Order,Type)                 \
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
typedef Cthermo<Mixture<icoPoly##Order##Type##ThermoPhysics> >                \
    Cthermo##Mixture##icoPoly##Order##Type##ThermoPhysics;                    \
                                                                              \
defineTemplateTypeNameAndDebugWithName                                        \
(                                                                             \
    Cthermo##Mixture##icoPoly##Order##Type##ThermoPhysics,                    \
    #Cthermo"<"#Mixture"<icoPoly"#Order#Type"ThermoPhysics>>",                \
    0                                                                         \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    BaseThermo,                                                               \
    Cthermo##Mixture##icoPoly##Order##Type##ThermoPhysics,                    \
    fvMesh                                                                    \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    fluidThermo,                                                              \
    Cthermo##Mixture##icoPoly##Order##Type##ThermoPhysics,                    \
    fvMesh                                                                    \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
