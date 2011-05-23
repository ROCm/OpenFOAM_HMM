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

#ifndef makeSprayParcelCollisionModels_H
#define makeSprayParcelCollisionModels_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoPhysicsTypes.H"
#include "SprayCloud.H"

#include "NoCollision.H"
#include "ORourkeCollision.H"
#include "TrajectoryCollision.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeSprayCollisionModels(ParcelType)                                  \
                                                                              \
    makeSprayCollisionModelThermoType                                         \
    (                                                                         \
        ParcelType,                                                           \
        constGasThermoPhysics                                                 \
    );                                                                        \
    makeSprayCollisionModelThermoType                                         \
    (                                                                         \
        ParcelType,                                                           \
        gasThermoPhysics                                                      \
    );                                                                        \
    makeSprayCollisionModelThermoType                                         \
    (                                                                         \
        ParcelType,                                                           \
        icoPoly8ThermoPhysics                                                 \
    );

#define makeSprayCollisionModelThermoType(ParcelType, ThermoType)             \
                                                                              \
    makeCollisionModel(SprayCloud<ParcelType<ThermoType> >);                  \
                                                                              \
    makeCollisionModelThermoType                                              \
    (                                                                         \
        NoCollision,                                                          \
        SprayCloud,                                                           \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
                                                                              \
    makeCollisionModelThermoType                                              \
    (                                                                         \
        ORourkeCollision,                                                     \
        SprayCloud,                                                           \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
                                                                              \
    makeCollisionModelThermoType                                              \
    (                                                                         \
        TrajectoryCollision,                                                  \
        SprayCloud,                                                           \
        ParcelType,                                                           \
        ThermoType                                                            \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
