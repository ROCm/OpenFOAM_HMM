/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#ifndef makeCoalParcelSurfaceReactionModels_H
#define makeCoalParcelSurfaceReactionModels_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "reactingThermoTypes.H"
#include "ReactingMultiphaseCloud.H"

#include "NoSurfaceReaction.H"
#include "COxidationDiffusionLimitedRate.H"
#include "COxidationKineticDiffusionLimitedRate.H"
#include "COxidationMurphyShaddix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define createCoalSurfaceReactionModelType(ParcelType)                        \
                                                                              \
    createCoalSurfaceReactionModelThermoType                                  \
    (                                                                         \
        ParcelType,                                                           \
        specieConstProperties                                                 \
    );                                                                        \
    createCoalSurfaceReactionModelThermoType                                  \
    (                                                                         \
        ParcelType,                                                           \
        specieReactingProperties                                              \
    );


#define createCoalSurfaceReactionModelThermoType(ParcelType, ThermoType)      \
                                                                              \
    makeSurfaceReactionModel                                                  \
    (                                                                         \
        ReactingMultiphaseCloud<ParcelType<ThermoType> >                      \
    );                                                                        \
                                                                              \
    makeSurfaceReactionModelThermoTypeNew                                     \
    (                                                                         \
        NoSurfaceReaction,                                                    \
        ReactingMultiphaseCloud,                                              \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
    makeSurfaceReactionModelThermoTypeNew                                     \
    (                                                                         \
        COxidationDiffusionLimitedRate,                                       \
        ReactingMultiphaseCloud,                                              \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
    makeSurfaceReactionModelThermoTypeNew                                     \
    (                                                                         \
        COxidationKineticDiffusionLimitedRate,                                \
        ReactingMultiphaseCloud,                                              \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
    makeSurfaceReactionModelThermoTypeNew                                     \
    (                                                                         \
        COxidationMurphyShaddix,                                              \
        ReactingMultiphaseCloud,                                              \
        ParcelType,                                                           \
        ThermoType                                                            \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
