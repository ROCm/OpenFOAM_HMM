/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "parcel.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// The diameter based Reynolds number
Foam::scalar Foam::parcel::Re
(
    const vector& U,
    const scalar nu
) const
{
    return mag(Urel(U))*d_/nu;
}

// The diameter based Reynolds number
Foam::scalar Foam::parcel::Re
(
    const scalar rho,
    const vector& U,
    const scalar mu
) const
{

    return rho*mag(Urel(U))*d_/mu;
}

// The diameter based Weber number
Foam::scalar Foam::parcel::We
(
    const vector& U,
    const scalar rho,
    const scalar sigma
) const
{
    return 0.5*rho*pow(mag(Urel(U)),2)*d_/sigma;
}


Foam::scalar Foam::parcel::Sc
(
    const scalar mu,
    const scalar rho,
    const scalar massDiffusion
) const
{
    return mu/(rho*massDiffusion);
}


Foam::scalar Foam::parcel::Sc
(
    const scalar nu,
    const scalar massDiffusion
) const
{
    return nu/massDiffusion;
}


Foam::scalar Foam::parcel::Pr
(
    const scalar cp,
    const scalar mu,
    const scalar kappa
) const
{
    return cp*mu/kappa;
}


Foam::scalar Foam::parcel::N(const scalar rho) const
{
    return 6.0*m_/(rho*pow3(d_)*constant::mathematical::pi);
}


Foam::scalar Foam::parcel::Vd() const
{
    return pow3(d_)*constant::mathematical::pi/6.0;
}


Foam::scalar Foam::parcel::V(const scalar rho) const
{
    return m_/rho;
}


// ************************************************************************* //
