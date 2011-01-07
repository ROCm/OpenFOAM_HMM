/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

#ifndef createKinematicParcelTypes_H
#define createKinematicParcelTypes_H

#include "KinematicParcel.H"
#include "KinematicCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define createKinematicParcelTypes(ParcelType)                                \
                                                                              \
    defineTypeNameAndDebug(ParcelType, 0);                                    \
    defineTemplateTypeNameAndDebug(Particle<ParcelType>, 0);                  \
    defineTemplateTypeNameAndDebug(Cloud<ParcelType>, 0);                     \
                                                                              \
    defineParcelTypeNameAndDebug(KinematicParcel<ParcelType>, 0);             \
    defineTemplateTypeNameAndDebug(KinematicParcel<ParcelType>, 0);           \
    defineParcelTypeNameAndDebug(KinematicCloud<ParcelType>, 0);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

