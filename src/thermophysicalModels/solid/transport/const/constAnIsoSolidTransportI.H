/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class thermo>
inline Foam::constAnIsoSolidTransport<thermo>::constAnIsoSolidTransport
(
    const thermo& t,
    const vector kappa
)
:
    thermo(t),
    kappa_(kappa)
{}


template<class thermo>
inline Foam::constAnIsoSolidTransport<thermo>::constAnIsoSolidTransport
(
    const word& name,
    const constAnIsoSolidTransport& ct
)
:
    thermo(name, ct),
    kappa_(ct.kappa_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class thermo>
inline Foam::scalar Foam::constAnIsoSolidTransport<thermo>::
kappa(const scalar T) const
{
    return mag(kappa_);
}

template<class thermo>
inline Foam::vector Foam::constAnIsoSolidTransport<thermo>::
Kappa(const scalar T) const
{
    return kappa_;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class thermo>
inline Foam::constAnIsoSolidTransport<thermo>&
Foam::constAnIsoSolidTransport<thermo>::operator=
(
    const constAnIsoSolidTransport<thermo>& ct
)
{
    kappa_ = ct.kappa_;

    return *this;
}


template<class thermo>
inline void Foam::constAnIsoSolidTransport<thermo>::operator+=
(
    const constAnIsoSolidTransport<thermo>& ct
)
{
    scalar molr1 = this->nMoles();

    molr1 /= this->nMoles();
    scalar molr2 = ct.nMoles()/this->nMoles();

    kappa_ = molr1*kappa_ + molr2*ct.kappa_;
}


template<class thermo>
inline void Foam::constAnIsoSolidTransport<thermo>::operator-=
(
    const constAnIsoSolidTransport<thermo>& ct
)
{
    scalar molr1 = this->nMoles();

    molr1 /= this->nMoles();
    scalar molr2 = ct.nMoles()/this->nMoles();

    kappa_ = molr1*kappa_ - molr2*ct.kappa_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class thermo>
inline Foam::constAnIsoSolidTransport<thermo> Foam::operator*
(
    const scalar s,
    const constAnIsoSolidTransport<thermo>& ct
)
{
    return constAnIsoSolidTransport<thermo>
    (
        s*static_cast<const thermo&>(ct),
        ct.kappa_
    );
}

// ************************************************************************* //
