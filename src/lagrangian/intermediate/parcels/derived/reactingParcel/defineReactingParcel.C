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

#include "reactingParcel.H"
#include "ReactingCloud.H"

namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<reactingParcel>, 0);

    defineParcelTypeNameAndDebug(KinematicParcel<reactingParcel>, 0);
//    defineTemplateTypeNameAndDebug(KinematicParcel<reactingParcel>, 0);
    defineParcelTypeNameAndDebug(ThermoParcel<reactingParcel>, 0);
    defineTemplateTypeNameAndDebug(ThermoParcel<reactingParcel>, 0);
    defineParcelTypeNameAndDebug(ReactingParcel<reactingParcel>, 0);
    defineTemplateTypeNameAndDebug(ReactingParcel<reactingParcel>, 0);

    defineParcelTypeNameAndDebug(KinematicCloud<reactingParcel>, 0);
//    defineTemplateTypeNameAndDebug(KinematicCloud<reactingParcel>, 0);

    defineParcelTypeNameAndDebug(ThermoCloud<reactingParcel>, 0);
//    defineTemplateTypeNameAndDebug(ThermoCloud<reactingParcel>, 0);

    defineParcelTypeNameAndDebug(ReactingCloud<reactingParcel>, 0);
//    defineTemplateTypeNameAndDebug(ReactingCloud<reactingParcel>, 0);

};


// ************************************************************************* //
