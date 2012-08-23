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

#include "constSolidRad.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class thermo>
inline Foam::constSolidRad<thermo>::constSolidRad
(
    const thermo& t,
    const scalar kappaRad,
    const scalar sigmaS,
    const scalar emissivity
)
:
    thermo(t),
    kappaRad_(kappaRad),
    sigmaS_(sigmaS),
    emissivity_(emissivity)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class thermo>
inline Foam::constSolidRad<thermo>::constSolidRad
(
    const word& name,
    const constSolidRad& pg
)
:
    thermo(name, pg),
    kappaRad_(pg.kappaRad_),
    sigmaS_(pg.sigmaS_),
    emissivity_(pg.emissivity_)
{}


template<class thermo>
inline Foam::autoPtr<Foam::constSolidRad<thermo> >
Foam::constSolidRad<thermo>::clone() const
{
    return autoPtr<constSolidRad<thermo> >(new constSolidRad<thermo>(*this));
}


template<class thermo>
inline Foam::autoPtr<Foam::constSolidRad<thermo> >
Foam::constSolidRad<thermo>::New(Istream& is)
{
    return autoPtr<constSolidRad<thermo> >(new constSolidRad<thermo>(is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class thermo>
inline Foam::scalar Foam::constSolidRad<thermo>::kappaRad(scalar) const
{
    return kappaRad_;
}


template<class thermo>
inline Foam::scalar Foam::constSolidRad<thermo>::sigmaS(scalar) const
{
    return sigmaS_;
}


template<class thermo>
inline Foam::scalar Foam::constSolidRad<thermo>::emissivity(scalar) const
{
    return emissivity_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class thermo>
inline Foam::constSolidRad<thermo>&
Foam::constSolidRad<thermo>::operator=(const constSolidRad<thermo>& ct)
{
    thermo::operator=(ct);

    kappaRad_ = ct.kappaRad_;
    sigmaS_ = ct.sigmaS_;
    emissivity_ = ct.emissivity_;

    return *this;
}


template<class thermo>
inline void Foam::constSolidRad<thermo>::operator+=
(
    const constSolidRad<thermo>& ct
)
{
    scalar molr1 = this->nMoles();

    thermo::operator+=(ct);

    molr1 /= this->nMoles();
    scalar molr2 = ct.nMoles()/this->nMoles();

    kappaRad_ = molr1*kappaRad_ + molr2*ct.kappaRad_;
    sigmaS_ = molr1*sigmaS_ + molr2*ct.sigmaS_;
    emissivity_ = molr1*emissivity_ + molr2*ct.emissivity_;
}


template<class thermo>
inline void Foam::constSolidRad<thermo>::operator-=
(
    const constSolidRad<thermo>& ct
)
{
    scalar molr1 = this->nMoles();

    thermo::operator-=(ct);

    molr1 /= this->nMoles();
    scalar molr2 = ct.nMoles()/this->nMoles();

    kappaRad_ = molr1*kappaRad_ - molr2*ct.kappaRad_;
    sigmaS_ = molr1*sigmaS_ - molr2*ct.sigmaS_;
    emissivity_ = molr1*emissivity_ - molr2*ct.emissivity_;
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class thermo>
inline Foam::constSolidRad<thermo> Foam::operator*
(
    const scalar s,
    const constSolidRad<thermo>& ct
)
{
    return constSolidRad<thermo>
    (
        s*static_cast<const thermo&>(ct),
        ct.kappaRad_,
        ct.sigmaS_,
        ct.emissivity_
    );
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// ************************************************************************* //
