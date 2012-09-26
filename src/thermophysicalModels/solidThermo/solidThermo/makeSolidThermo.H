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

#define makeSolidThermo(BaseThermo,Cthermo,Mixture,Transport,Radiation,Type,Thermo,EqnOfState,Specie)\
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
                        EqnOfState                                            \
                        <                                                     \
                            Specie                                            \
                        >                                                     \
                    >,                                                        \
                    Type                                                      \
                >                                                             \
            >                                                                 \
        >                                                                     \
    >,                                                                        \
    BaseThermo                                                                \
> Cthermo##Mixture##Transport##Radiation##Type##Thermo                        \
      ##EqnOfState##Specie##BaseThermo;                                       \
                                                                              \
defineTemplateTypeNameAndDebugWithName                                        \
(                                                                             \
    Cthermo##Mixture##Transport##Radiation##Type##Thermo                      \
         ##EqnOfState##Specie##BaseThermo,                                    \
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
                        #EqnOfState                                           \
                        "<"                                                   \
                            #Specie                                           \
                        ">"                                                   \
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
                    EqnOfState                                                \
                    <                                                         \
                        Specie                                                \
                    >                                                         \
                >,                                                            \
                Type                                                          \
            >                                                                 \
        >                                                                     \
    >                                                                         \
> Mixture##Transport##Radiation##Type##Thermo##EqnOfState##Specie;            \
                                                                              \
                                                                              \
defineTemplateTypeNameAndDebugWithName                                        \
(                                                                             \
    Mixture##Transport##Radiation##Type##Thermo##EqnOfState##Specie,          \
    #Mixture                                                                  \
    "<"                                                                       \
        #Transport                                                            \
        "<"                                                                   \
            #Radiation                                                        \
            "<"                                                               \
                #Thermo                                                       \
                "<"                                                           \
                    #EqnOfState                                               \
                ">,"                                                          \
                #Type                                                         \
    ">>>",                                                                    \
    0                                                                         \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    BaseThermo,                                                               \
    Cthermo##Mixture##Transport##Radiation##Type##Thermo                      \
        ##EqnOfState##Specie##BaseThermo,                                     \
    mesh                                                                      \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    basicThermo,                                                              \
    Cthermo##Mixture##Transport##Radiation##Type##Thermo                      \
        ##EqnOfState##Specie##BaseThermo,                                     \
    fvMesh                                                                    \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    BaseThermo,                                                               \
    Cthermo##Mixture##Transport##Radiation##Type##Thermo                      \
        ##EqnOfState##Specie##BaseThermo,                                     \
    dictionary                                                                \
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
