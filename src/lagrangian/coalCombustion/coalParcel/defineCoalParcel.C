/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "coalParcel.H"
#include "ReactingMultiphaseCloud.H"

namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<coalParcel>, 0);

    defineParcelTypeNameAndDebug(KinematicParcel<coalParcel>, 0);
//    defineTemplateTypeNameAndDebug(KinematicParcel<coalParcel>, 0);
    defineParcelTypeNameAndDebug(ThermoParcel<coalParcel>, 0);
    defineTemplateTypeNameAndDebug(ThermoParcel<coalParcel>, 0);
    defineParcelTypeNameAndDebug(ReactingParcel<coalParcel>, 0);
    defineTemplateTypeNameAndDebug(ReactingParcel<coalParcel>, 0);
    defineParcelTypeNameAndDebug(ReactingMultiphaseParcel<coalParcel>, 0);
    defineTemplateTypeNameAndDebug(ReactingMultiphaseParcel<coalParcel>, 0);

    defineParcelTypeNameAndDebug(KinematicCloud<coalParcel>, 0);
    defineParcelTypeNameAndDebug(ThermoCloud<coalParcel>, 0);
    defineParcelTypeNameAndDebug(ReactingCloud<coalParcel>, 0);
    defineParcelTypeNameAndDebug(ReactingMultiphaseCloud<coalParcel>, 0);
};


// ************************************************************************* //
