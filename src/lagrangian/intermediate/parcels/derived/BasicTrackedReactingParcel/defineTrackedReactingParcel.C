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

#include "trackedReactingParcel.H"
#include "ReactingCloud.H"

namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<trackedReactingParcel>, 0);

    defineParcelTypeNameAndDebug(KinematicParcel<trackedReactingParcel>, 0);
//    defineTemplateTypeNameAndDebug(KinematicParcel<trackedReactingParcel>, 0);
    defineParcelTypeNameAndDebug(ThermoParcel<trackedReactingParcel>, 0);
    defineTemplateTypeNameAndDebug(ThermoParcel<trackedReactingParcel>, 0);
    defineParcelTypeNameAndDebug(ReactingParcel<trackedReactingParcel>, 0);
    defineTemplateTypeNameAndDebug(ReactingParcel<trackedReactingParcel>, 0);
    defineTemplateTypeNameAndDebug
    (
        TrackedReactingParcel<trackedReactingParcel>,
        0
    );

    defineParcelTypeNameAndDebug(KinematicCloud<trackedReactingParcel>, 0);
//    defineTemplateTypeNameAndDebug(KinematicCloud<trackedReactingParcel>, 0);

    defineParcelTypeNameAndDebug(ThermoCloud<trackedReactingParcel>, 0);
//    defineTemplateTypeNameAndDebug(ThermoCloud<trackedReactingParcel>, 0);

    defineParcelTypeNameAndDebug(ReactingCloud<trackedReactingParcel>, 0);
//    defineTemplateTypeNameAndDebug(ReactingCloud<trackedReactingParcel>, 0);
}


// ************************************************************************* //
