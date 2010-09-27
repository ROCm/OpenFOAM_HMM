/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "basicReactingParcel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicReactingParcel::basicReactingParcel
(
    ReactingCloud<basicReactingParcel>& owner,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI
)
:
    ReactingParcel<basicReactingParcel>
    (
        owner,
        position,
        cellI,
        tetFaceI,
        tetPtI
    )
{}


Foam::basicReactingParcel::basicReactingParcel
(
    ReactingCloud<basicReactingParcel>& owner,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const scalar dTarget0,
    const vector& U0,
    const vector& f0,
    const vector& angularMomentum0,
    const vector& torque0,
    const scalarField& Y0,
    const ReactingParcel<basicReactingParcel>::constantProperties& constProps
)
:
    ReactingParcel<basicReactingParcel>
    (
        owner,
        position,
        cellI,
        tetFaceI,
        tetPtI,
        typeId,
        nParticle0,
        d0,
        dTarget0,
        U0,
        f0,
        angularMomentum0,
        torque0,
        Y0,
        constProps
    )
{}


Foam::basicReactingParcel::basicReactingParcel
(
    const Cloud<basicReactingParcel>& cloud,
    Istream& is,
    bool readFields
)
:
    ReactingParcel<basicReactingParcel>(cloud, is, readFields)
{}


Foam::basicReactingParcel::basicReactingParcel
(
    const basicReactingParcel& p
)
:
    ReactingParcel<basicReactingParcel>(p)
{}


// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //

Foam::basicReactingParcel::~basicReactingParcel()
{}


// ************************************************************************* //
