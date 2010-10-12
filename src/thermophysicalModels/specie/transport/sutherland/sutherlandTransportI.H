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

#include "specie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class thermo>
inline void sutherlandTransport<thermo>::calcCoeffs
(
    const scalar mu1, const scalar T1,
    const scalar mu2, const scalar T2
)
{
    scalar rootT1 = sqrt(T1);
    scalar mu1rootT2 = mu1*sqrt(T2);
    scalar mu2rootT1 = mu2*rootT1;

    Ts_ = (mu2rootT1 - mu1rootT2)/(mu1rootT2/T1 - mu2rootT1/T2);

    As_ = mu1*(1.0 + Ts_/T1)/rootT1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class thermo>
inline sutherlandTransport<thermo>::sutherlandTransport
(
    const thermo& t,
    const scalar as,
    const scalar ts
)
:
    thermo(t),
    As_(as),
    Ts_(ts)
{}


template<class thermo>
inline sutherlandTransport<thermo>::sutherlandTransport
(
    const thermo& t,
    const scalar mu1, const scalar T1,
    const scalar mu2, const scalar T2
)
:
    thermo(t)
{
    calcCoeffs(mu1, T1, mu2, T2);
}


template<class thermo>
inline sutherlandTransport<thermo>::sutherlandTransport
(
    const word& name,
    const sutherlandTransport& st
)
:
    thermo(name, st),
    As_(st.As_),
    Ts_(st.Ts_)
{}


template<class thermo>
inline autoPtr<sutherlandTransport<thermo> > sutherlandTransport<thermo>::clone
() const
{
    return autoPtr<sutherlandTransport<thermo> >
    (
        new sutherlandTransport<thermo>(*this)
    );
}


template<class thermo>
inline autoPtr<sutherlandTransport<thermo> > sutherlandTransport<thermo>::New
(
    Istream& is
)
{
    return autoPtr<sutherlandTransport<thermo> >
    (
        new sutherlandTransport<thermo>(is)
    );
}


template<class thermo>
inline autoPtr<sutherlandTransport<thermo> > sutherlandTransport<thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<sutherlandTransport<thermo> >
    (
        new sutherlandTransport<thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Dynamic viscosity [kg/ms]
template<class thermo>
inline scalar sutherlandTransport<thermo>::mu(const scalar T) const
{
    return As_*::sqrt(T)/(1.0 + Ts_/T);
}


template<class thermo>
inline scalar sutherlandTransport<thermo>::kappa(const scalar T) const
{
    scalar Cv_ = this->Cv(T);
    return mu(T)*Cv_*(1.32 + 1.77*this->R()/Cv_);
}


template<class thermo>
inline scalar sutherlandTransport<thermo>::alpha(const scalar T) const
{
    scalar Cv_ = this->Cv(T);
    scalar R_ = this->R();
    scalar Cp_ = Cv_ + R_;

    scalar deltaT = T - specie::Tstd;
    scalar CpBar =
        (deltaT*(this->H(T) - this->H(specie::Tstd)) + Cp_)/(sqr(deltaT) + 1);

    return mu(T)*Cv_*(1.32 + 1.77*this->R()/Cv_)/CpBar;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class thermo>
inline sutherlandTransport<thermo>& sutherlandTransport<thermo>::operator=
(
    const sutherlandTransport<thermo>& st
)
{
    thermo::operator=(st);

    As_ = st.As_;
    Ts_ = st.Ts_;

    return *this;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class thermo>
inline sutherlandTransport<thermo> operator+
(
    const sutherlandTransport<thermo>& st1,
    const sutherlandTransport<thermo>& st2
)
{
    thermo t
    (
        static_cast<const thermo&>(st1) + static_cast<const thermo&>(st2)
    );

    scalar molr1 = st1.nMoles()/t.nMoles();
    scalar molr2 = st2.nMoles()/t.nMoles();

    return sutherlandTransport<thermo>
    (
        t,
        molr1*st1.As_ + molr2*st2.As_,
        molr1*st1.Ts_ + molr2*st2.Ts_
    );
}


template<class thermo>
inline sutherlandTransport<thermo> operator-
(
    const sutherlandTransport<thermo>& st1,
    const sutherlandTransport<thermo>& st2
)
{
    thermo t
    (
        static_cast<const thermo&>(st1) - static_cast<const thermo&>(st2)
    );

    scalar molr1 = st1.nMoles()/t.nMoles();
    scalar molr2 = st2.nMoles()/t.nMoles();

    return sutherlandTransport<thermo>
    (
        t,
        molr1*st1.As_ - molr2*st2.As_,
        molr1*st1.Ts_ - molr2*st2.Ts_
    );
}


template<class thermo>
inline sutherlandTransport<thermo> operator*
(
    const scalar s,
    const sutherlandTransport<thermo>& st
)
{
    return sutherlandTransport<thermo>
    (
        s*static_cast<const thermo&>(st),
        st.As_,
        st.Ts_
    );
}


template<class thermo>
inline sutherlandTransport<thermo> operator==
(
    const sutherlandTransport<thermo>& st1,
    const sutherlandTransport<thermo>& st2
)
{
    return st2 - st1;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
