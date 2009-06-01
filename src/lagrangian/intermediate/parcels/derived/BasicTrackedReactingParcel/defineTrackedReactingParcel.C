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

#include "basicTrackedReactingParcelTypes.H"
#include "BasicReactingCloud.H"

namespace Foam
{
    defineTemplateTypeNameAndDebug(bTrackedReactingParcel, 0);
    defineTemplateTypeNameAndDebug(Particle<bTrackedReactingParcel>, 0);
    defineTemplateTypeNameAndDebug(Cloud<bTrackedReactingParcel>, 0);

    defineParcelTypeNameAndDebug(KinematicParcel<bTrackedReactingParcel>, 0);
    defineTemplateTypeNameAndDebug(KinematicParcel<bTrackedReactingParcel>, 0);
    defineParcelTypeNameAndDebug(ThermoParcel<bTrackedReactingParcel>, 0);
    defineTemplateTypeNameAndDebug(ThermoParcel<bTrackedReactingParcel>, 0);
    defineParcelTypeNameAndDebug(ReactingParcel<bTrackedReactingParcel>, 0);
    defineTemplateTypeNameAndDebug(ReactingParcel<bTrackedReactingParcel>, 0);
    defineTemplateTypeNameAndDebug
    (
        TrackedReactingParcel<bTrackedReactingParcel>,
        0
    );

    defineParcelTypeNameAndDebug(KinematicCloud<bTrackedReactingParcel>, 0);
//    defineTemplateTypeNameAndDebug(KinematicCloud<bTrackedReactingParcel>, 0);

    defineParcelTypeNameAndDebug(ThermoCloud<bTrackedReactingParcel>, 0);
//    defineTemplateTypeNameAndDebug(ThermoCloud<bTrackedReactingParcel>, 0);

    defineParcelTypeNameAndDebug(ReactingCloud<bTrackedReactingParcel>, 0);
//    defineTemplateTypeNameAndDebug(ReactingCloud<bTrackedReactingParcel>, 0);
}


// ************************************************************************* //
