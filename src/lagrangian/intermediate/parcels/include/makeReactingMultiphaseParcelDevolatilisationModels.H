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

#ifndef makeReactingMultiphaseParcelDevolatilisationModels_H
#define makeReactingMultiphaseParcelDevolatilisationModels_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Cloud.H"
#include "KinematicCloud.H"
#include "ThermoCloud.H"
#include "ReactingCloud.H"
#include "ReactingMultiphaseCloud.H"

#include "ConstantRateDevolatilisation.H"
#include "NoDevolatilisation.H"
#include "SingleKineticRateDevolatilisation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeReactingMultiphaseParcelDevolatilisationModels(ParcelType)        \
                                                                              \
    typedef ReactingMultiphaseCloud                                           \
    <                                                                         \
        ReactingCloud                                                         \
        <                                                                     \
            ThermoCloud                                                       \
            <                                                                 \
                KinematicCloud                                                \
                <                                                             \
                    Cloud                                                     \
                    <                                                         \
                        ParcelType                                            \
                    >                                                         \
                >                                                             \
            >                                                                 \
        >                                                                     \
    > RMPRTKC##ParcelType;                                                    \
                                                                              \
    makeDevolatilisationModel(RMPRTKC##ParcelType);                           \
    makeDevolatilisationModelType                                             \
    (                                                                         \
        ConstantRateDevolatilisation,                                         \
        RMPRTKC##ParcelType                                                   \
    );                                                                        \
    makeDevolatilisationModelType                                             \
    (                                                                         \
        NoDevolatilisation,                                                   \
        RMPRTKC##ParcelType                                                   \
    );                                                                        \
    makeDevolatilisationModelType                                             \
    (                                                                         \
        SingleKineticRateDevolatilisation,                                    \
        RMPRTKC##ParcelType                                                   \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
