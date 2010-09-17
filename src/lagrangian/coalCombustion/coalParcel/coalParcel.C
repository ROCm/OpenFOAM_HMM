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

#include "coalParcel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coalParcel::coalParcel
(
    ReactingMultiphaseCloud<coalParcel>& owner,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI
)
:
    ReactingMultiphaseParcel<coalParcel>
    (
        owner,
        position,
        cellI,
        tetFaceI,
        tetPtI
    )
{}


Foam::coalParcel::coalParcel
(
    ReactingMultiphaseCloud<coalParcel>& owner,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const vector& U0,
    const vector& f0,
    const vector& pi0,
    const vector& tau0,
    const scalarField& YMixture0,
    const scalarField& YGas0,
    const scalarField& YLiquid0,
    const scalarField& YSolid0,
    const ReactingMultiphaseParcel<coalParcel>::constantProperties& constProps
)
:
    ReactingMultiphaseParcel<coalParcel>
    (
        owner,
        position,
        cellI,
        tetFaceI,
        tetPtI,
        typeId,
        nParticle0,
        d0,
        U0,
        f0,
        pi0,
        tau0,
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


Foam::coalParcel::coalParcel(const coalParcel& p)
:
    ReactingMultiphaseParcel<coalParcel>(p)
{
}


// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //

Foam::coalParcel::~coalParcel()
{}


// ************************************************************************* //
