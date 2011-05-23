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

#ifndef makeSprayParcelBreakupModels_H
#define makeSprayParcelBreakupModels_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoPhysicsTypes.H"
#include "SprayCloud.H"

#include "NoBreakup.H"
#include "ReitzDiwakar.H"
#include "ReitzKHRT.H"
#include "TAB.H"
#include "ETAB.H"
#include "SHF.H"
#include "PilchErdman.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeSprayBreakupModels(ParcelType)                                    \
                                                                              \
    makeSprayBreakupModelThermoType                                           \
    (                                                                         \
        ParcelType,                                                           \
        constGasThermoPhysics                                                 \
    );                                                                        \
                                                                              \
    makeSprayBreakupModelThermoType                                           \
    (                                                                         \
        ParcelType,                                                           \
        gasThermoPhysics                                                      \
    );                                                                        \
                                                                              \
    makeSprayBreakupModelThermoType                                           \
    (                                                                         \
        ParcelType,                                                           \
        icoPoly8ThermoPhysics                                                 \
    );

#define makeSprayBreakupModelThermoType(ParcelType, ThermoType)               \
                                                                              \
    makeBreakupModel(SprayCloud<ParcelType<ThermoType> >);                    \
                                                                              \
    makeBreakupModelThermoType                                                \
    (                                                                         \
        NoBreakup,                                                            \
        SprayCloud,                                                           \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
                                                                              \
    makeBreakupModelThermoType                                                \
    (                                                                         \
        PilchErdman,                                                          \
        SprayCloud,                                                           \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
                                                                              \
    makeBreakupModelThermoType                                                \
    (                                                                         \
        SHF,                                                                  \
        SprayCloud,                                                           \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
                                                                              \
    makeBreakupModelThermoType                                                \
    (                                                                         \
        ETAB,                                                                 \
        SprayCloud,                                                           \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
                                                                              \
    makeBreakupModelThermoType                                                \
    (                                                                         \
        TAB,                                                                  \
        SprayCloud,                                                           \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
                                                                              \
    makeBreakupModelThermoType                                                \
    (                                                                         \
        ReitzKHRT,                                                            \
        SprayCloud,                                                           \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
                                                                              \
    makeBreakupModelThermoType                                                \
    (                                                                         \
        ReitzDiwakar,                                                         \
        SprayCloud,                                                           \
        ParcelType,                                                           \
        ThermoType                                                            \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
