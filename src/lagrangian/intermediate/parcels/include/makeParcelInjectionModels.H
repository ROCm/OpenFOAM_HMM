/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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

#ifndef makeParcelInjectionModels_H
#define makeParcelInjectionModels_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "InteractingKinematicCloud.H"

#include "ConeInjection.H"
#include "ConeInjectionMP.H"
#include "FieldActivatedInjection.H"
#include "KinematicLookupTableInjection.H"
#include "ManualInjection.H"
#include "NoInjection.H"
#include "PatchInjection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeParcelInjectionModels(ParcelType)                                 \
                                                                              \
    makeInjectionModel(InteractingKinematicCloud<ParcelType>);                           \
                                                                              \
    makeInjectionModelType                                                    \
    (                                                                         \
        ConeInjection,                                                        \
        InteractingKinematicCloud,                                                       \
        ParcelType                                                            \
    );                                                                        \
    makeInjectionModelType                                                    \
    (                                                                         \
        ConeInjectionMP,                                                      \
        InteractingKinematicCloud,                                                       \
        ParcelType                                                            \
    );                                                                        \
    makeInjectionModelType                                                    \
    (                                                                         \
        FieldActivatedInjection,                                              \
        InteractingKinematicCloud,                                                       \
        ParcelType                                                            \
    );                                                                        \
    makeInjectionModelType                                                    \
    (                                                                         \
        KinematicLookupTableInjection,                                        \
        InteractingKinematicCloud,                                                       \
        ParcelType                                                            \
    );                                                                        \
    makeInjectionModelType                                                    \
    (                                                                         \
        ManualInjection,                                                      \
        InteractingKinematicCloud,                                                       \
        ParcelType                                                            \
    );                                                                        \
    makeInjectionModelType                                                    \
    (                                                                         \
        NoInjection,                                                          \
        InteractingKinematicCloud,                                                       \
        ParcelType                                                            \
    );                                                                        \
    makeInjectionModelType                                                    \
    (                                                                         \
        PatchInjection,                                                       \
        InteractingKinematicCloud,                                                       \
        ParcelType                                                            \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
