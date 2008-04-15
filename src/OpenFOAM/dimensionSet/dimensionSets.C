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

Description
    Useful dimension sets.

\*---------------------------------------------------------------------------*/

#include "dimensionSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const dimensionSet dimless(0, 0, 0, 0, 0, 0, 0);

const dimensionSet dimMass(1, 0, 0, 0, 0, 0, 0);
const dimensionSet dimLength(0, 1, 0, 0, 0, 0, 0);
const dimensionSet dimTime(0, 0, 1, 0, 0, 0, 0);
const dimensionSet dimTemperature(0, 0, 0, 1, 0, 0, 0);
const dimensionSet dimMoles(0, 0, 0, 0, 1, 0, 0);
const dimensionSet dimCurrent(0, 0, 0, 0, 0, 1, 0);
const dimensionSet dimLuminousIntensity(0, 0, 0, 0, 0, 0, 1);

const dimensionSet dimArea(0, 2, 0, 0, 0, 0, 0);
const dimensionSet dimVolume(0, 3, 0, 0, 0, 0, 0);
const dimensionSet dimVol(0, 3, 0, 0, 0, 0, 0);

const dimensionSet dimDensity(1, -3, 0, 0, 0, 0, 0);
const dimensionSet dimForce(1, 1, -2, 0, 0, 0, 0);

const dimensionSet dimVelocity(0, 1, -1, 0, 0, 0, 0);
const dimensionSet dimPressure(1, -1, -2, 0, 0, 0, 0);
const dimensionSet dimGasConstant(0, 2, -2, -1, 0, 0, 0);
const dimensionSet dimSpecificHeatCapacity(0, 2, -2, -1, 0, 0, 0);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
