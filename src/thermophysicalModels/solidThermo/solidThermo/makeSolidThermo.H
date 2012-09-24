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
    Foam::solidThermo

Description
    Macros for creating solid thermo packages

\*---------------------------------------------------------------------------*/

#ifndef makeSolidThermo_H
#define makesolidThermo_H

#include "addToRunTimeSelectionTable.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

#define makeSolidThermo(BaseThermo,Cthermo,Mixture,Transport,Radiation,Type,Thermo,Rho)\
                                                                              \
typedef Cthermo                                                               \
<                                                                             \
    Mixture                                                                   \
    <                                                                         \
        Transport                                                             \
        <                                                                     \
            Radiation                                                         \
            <                                                                 \
                species::thermo                                               \
                <                                                             \
                    Thermo                                                    \
                    <                                                         \
                        Rho                                                   \
                    >,                                                        \
                    Type                                                      \
                >                                                             \
            >                                                                 \
        >                                                                     \
    >,                                                                        \
    BaseThermo                                                                \
>   Cthermo##Mixture##Transport##Radiation##Type##Thermo##Rho##BaseThermo;    \
                                                                              \
defineTemplateTypeNameAndDebugWithName                                        \
(                                                                             \
    Cthermo##Mixture##Transport##Radiation##Type##Thermo##Rho##BaseThermo,    \
    #Cthermo                                                                  \
    "<"                                                                       \
        #Mixture                                                              \
        "<"                                                                   \
            #Transport                                                        \
            "<"                                                               \
                #Radiation                                                    \
                "<"                                                           \
                    #Thermo                                                   \
                    "<"                                                       \
                        #Rho                                                  \
                    ">,"                                                      \
                    #Type                                                     \
        ">>>>",                                                               \
    0                                                                         \
);                                                                            \
                                                                              \
                                                                              \
typedef Mixture                                                               \
<                                                                             \
    Transport                                                                 \
    <                                                                         \
        Radiation                                                             \
        <                                                                     \
            species::thermo                                                   \
            <                                                                 \
                Thermo                                                        \
                <                                                             \
                    Rho                                                       \
                >,                                                            \
                Type                                                          \
            >                                                                 \
        >                                                                     \
    >                                                                         \
>   Mixture##Transport##Radiation##Type##Thermo##Rho;                         \
                                                                              \
                                                                              \
                                                                              \
defineTemplateTypeNameAndDebugWithName                                        \
(                                                                             \
    Mixture##Transport##Radiation##Type##Thermo##Rho,                         \
    #Mixture                                                                  \
    "<"                                                                       \
        #Transport                                                            \
        "<"                                                                   \
            #Radiation                                                        \
            "<"                                                               \
                #Thermo                                                       \
                "<"                                                           \
                    #Rho                                                      \
                ">,"                                                          \
                #Type                                                         \
    ">>>",                                                                    \
    0                                                                         \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    BaseThermo,                                                               \
    Cthermo##Mixture##Transport##Radiation##Type##Thermo##Rho##BaseThermo,    \
    mesh                                                                      \
);                                                                            \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    basicThermo,                                                              \
    Cthermo##Mixture##Transport##Radiation##Type##Thermo##Rho##BaseThermo,    \
    fvMesh                                                                    \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    BaseThermo,                                                               \
    Cthermo##Mixture##Transport##Radiation##Type##Thermo##Rho##BaseThermo,    \
    dictionary                                                                \
);


}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
