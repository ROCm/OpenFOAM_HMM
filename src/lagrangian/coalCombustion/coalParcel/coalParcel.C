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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coalParcel, 0);
    defineParticleTypeNameAndDebug(coalParcel, 0);
    defineParcelTypeNameAndDebug(coalParcel, 0);
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalParcel::coalParcel
(
    ReactingMultiphaseCloud<coalParcel>& owner,
    const vector& position,
    const label celli,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const vector& U0,
    const scalarField& YMixture0,
    const scalarField& YGas0,
    const scalarField& YLiquid0,
    const scalarField& YSolid0,
    const constantProperties& constProps
)
:
    ReactingMultiphaseParcel<coalParcel>
    (
        owner,
        position,
        celli,
        typeId,
        nParticle0,
        d0,
        U0,
        YMixture0,
        YGas0,
        YLiquid0,
        YSolid0,
        constProps
    )
{}


Foam::coalParcel::coalParcel
(
    const Cloud<coalParcel>& cloud,
    Istream& is,
    bool readFields
)
:
    ReactingMultiphaseParcel<coalParcel>(cloud, is, readFields)
{}


// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //

Foam::coalParcel::~coalParcel()
{}


// ************************************************************************* //
