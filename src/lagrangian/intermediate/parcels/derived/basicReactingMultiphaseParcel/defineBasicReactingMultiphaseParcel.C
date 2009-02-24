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

#include "basicReactingMultiphaseParcel.H"
#include "ReactingCloud.H"

namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<basicReactingMultiphaseParcel>, 0);

    defineParcelTypeNameAndDebug
    (
        KinematicParcel<basicReactingMultiphaseParcel>,
        0
    );
//    defineTemplateTypeNameAndDebug
//    (
//        KinematicParcel<basicReactingMultiphaseParcel>,
//        0
//    );
    defineParcelTypeNameAndDebug
    (
        ThermoParcel<basicReactingMultiphaseParcel>,
        0
    );
    defineTemplateTypeNameAndDebug
    (
        ThermoParcel<basicReactingMultiphaseParcel>,
        0
    );
    defineParcelTypeNameAndDebug
    (
        ReactingParcel<basicReactingMultiphaseParcel>,
        0
    );
    defineTemplateTypeNameAndDebug
    (
        ReactingParcel<basicReactingMultiphaseParcel>,
        0
    );

    defineParcelTypeNameAndDebug
    (
        KinematicCloud<basicReactingMultiphaseParcel>,
        0
    );
//    defineTemplateTypeNameAndDebug
//    (
//      KinematicCloud<basicReactingMultiphaseParcel>,
//      0
//    );

    defineParcelTypeNameAndDebug
    (
        ThermoCloud<basicReactingMultiphaseParcel>,
        0
    );
//    defineTemplateTypeNameAndDebug
//    (
//        ThermoCloud<basicReactingMultiphaseParcel>,
//        0
//    );

    defineParcelTypeNameAndDebug
    (
        ReactingCloud<basicReactingMultiphaseParcel>,
        0
    );
//    defineTemplateTypeNameAndDebug
//    (
//        ReactingCloud<basicReactingMultiphaseParcel>,
//        0
//    );

};


// ************************************************************************* //
