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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(trackedReactingParcel, 0);
    defineParticleTypeNameAndDebug(trackedReactingParcel, 0);
    defineParcelTypeNameAndDebug(trackedReactingParcel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::trackedReactingParcel::trackedReactingParcel
(
    ReactingCloud<trackedReactingParcel>& owner,
    const vector& position,
    const label cellI
)
:
    TrackedReactingParcel<trackedReactingParcel>(owner, position, cellI)
{}


Foam::trackedReactingParcel::trackedReactingParcel
(
    ReactingCloud<trackedReactingParcel>& owner,
    const vector& position,
    const label cellI,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const vector& U0,
    const scalarField& Y0,
    const constantProperties& constProps
)
:
    TrackedReactingParcel<trackedReactingParcel>
    (
        owner,
        position,
        cellI,
        typeId,
        nParticle0,
        d0,
        U0,
        Y0,
        constProps
    )
{}


Foam::trackedReactingParcel::trackedReactingParcel
(
    const Cloud<trackedReactingParcel>& cloud,
    Istream& is,
    bool readFields
)
:
    TrackedReactingParcel<trackedReactingParcel>(cloud, is, readFields)
{}


Foam::trackedReactingParcel::trackedReactingParcel
(
    const trackedReactingParcel& p
)
:
    TrackedReactingParcel<trackedReactingParcel>(p)
{}


// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //

Foam::trackedReactingParcel::~trackedReactingParcel()
{}


// ************************************************************************* //
