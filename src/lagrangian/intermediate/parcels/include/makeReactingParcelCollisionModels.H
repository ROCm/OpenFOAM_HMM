/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#ifndef makeReactingParcelCollisionModels_H
#define makeReactingParcelCollisionModels_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoPhysicsTypes.H"
#include "KinematicCloud.H"

#include "NoCollision.H"
#include "PairCollision.H"

#include "PairSpringSliderDashpot.H"

#include "WallSpringSliderDashpot.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeReactingCollisionModels(ParcelType)                               \
                                                                              \
    makeReactingCollisionModelThermoType                                      \
    (                                                                         \
        ParcelType,                                                           \
        constGasThermoPhysics                                                 \
    );                                                                        \
                                                                              \
    makeReactingCollisionModelThermoType                                      \
    (                                                                         \
        ParcelType,                                                           \
        gasThermoPhysics                                                      \
    );                                                                        \
                                                                              \
    makeReactingCollisionModelThermoType                                      \
    (                                                                         \
        ParcelType,                                                           \
        icoPoly8ThermoPhysics                                                 \
    );


#define makeReactingCollisionModelThermoType(ParcelType, ThermoType)          \
                                                                              \
    makeCollisionModel(KinematicCloud<ParcelType<ThermoType> >);              \
                                                                              \
    makeCollisionModelThermoType                                              \
    (                                                                         \
        NoCollision,                                                          \
        KinematicCloud,                                                       \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
                                                                              \
    makeCollisionModelThermoType                                              \
    (                                                                         \
        PairCollision,                                                        \
        KinematicCloud,                                                       \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
                                                                              \
    makePairModel(KinematicCloud<ParcelType<ThermoType> >);                   \
                                                                              \
    makePairModelThermoType                                                   \
    (                                                                         \
        PairSpringSliderDashpot,                                              \
        KinematicCloud,                                                       \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
                                                                              \
    makeWallModel(KinematicCloud<ParcelType<ThermoType> >);                   \
                                                                              \
    makeWallModelThermoType                                                   \
    (                                                                         \
        WallSpringSliderDashpot,                                              \
        KinematicCloud,                                                       \
        ParcelType,                                                           \
        ThermoType                                                            \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
