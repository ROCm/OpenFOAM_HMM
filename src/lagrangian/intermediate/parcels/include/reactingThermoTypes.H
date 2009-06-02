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

Typedefs
    Foam::ReactingCloudThermoTypes

Description
    Type definitions for available thermo pacakges

\*---------------------------------------------------------------------------*/

#ifndef reactingThermoTypes_H
#define reactingThermoTypes_H

#include "perfectGas.H"
#include "hConstThermo.H"
#include "janafThermo.H"
#include "specieThermo.H"
#include "sutherlandTransport.H"
#include "constTransport.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef sutherlandTransport<specieThermo<janafThermo<perfectGas> > >
        specieReactingProperties;

    typedef constTransport<specieThermo<hConstThermo<perfectGas> > >
        specieConstProperties;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
