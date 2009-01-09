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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::eConstThermo::eConstThermo
(
    const specieThermo& st,
    const scalar cv,
    const scalar hf
)
:
    specieThermo(st),
    Cv_(cv),
    Hf_(hf)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::eConstThermo::eConstThermo
(
    const word& name,
    const eConstThermo& ct
)
:
    specieThermo(name, ct),
    Cv_(ct.Cv_),
    Hf_(ct.Hf_)
{}


template<class equationOfState>
inline Foam::autoPtr<eConstThermo<equationOfState> >
Foam::eConstThermo<equationOfState>::clone() const
{
    return autoPtr<eConstThermo<equationOfState> >
    (
        new eConstThermo<equationOfState>(*this)
    );
}


template<class equationOfState>
inline autoPtr<eConstThermo<equationOfState> >
Foam::eConstThermo<equationOfState>::New(Istream& is)
{
    return autoPtr<eConstThermo<equationOfState> >
    (
        new eConstThermo<equationOfState>(is)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::eConstThermo::cp(const scalar) const
{
    return Cv_*W() + RR;
}


inline Foam::scalar Foam::eConstThermo::h(const scalar T) const
{
    return cp(T)*T + Hf_*W();
}


inline Foam::scalar Foam::eConstThermo::hs(const scalar T) const
{
    return cp(T)*T;
}


template<class equationOfState>
inline Foam::scalar Foam::eConstThermo<equationOfState>::hc() const
{
    return Hf_*this->W();
}


inline Foam::scalar Foam::eConstThermo::s(const scalar T) const
{
    notImplemented("scalar eConstThermo::s(const scalar T) const");
    return T;
}


inline Foam::scalar Foam::eConstThermo::TH
(
    const scalar h,
    const scalar T0
) const
{
    return (h - Hf_)/Cp(T0);
}


inline Foam::scalar Foam::eConstThermo::TE
(
    const scalar e,
    const scalar
) const
{
    return (e - Hf_)/Cv_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline Foam::eConstThermo& Foam::eConstThermo::operator=
(
    const eConstThermo& ct
)
{
    specieThermo::operator=(ct);

    Cv_ = ct.Cv_;
    Hf_ = ct.Hf_;

    return *this;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline Foam::eConstThermo Foam::operator+
(
    const eConstThermo& ct1,
    const eConstThermo& ct2
)
{
    specieThermo st(((const specieThermo&)ct1) + ((const specieThermo&)ct2));

    return eConstThermo
    (
        st,
        ct1.nMoles()/st.nMoles()*ct1.Cv_ + ct2.nMoles()/st.nMoles()*ct2.Cv_,
        ct1.nMoles()/st.nMoles()*ct1.Hf_ + ct2.nMoles()/st.nMoles()*ct2.Hf_
    );
}


inline Foam::eConstThermo Foam::operator-
(
    const eConstThermo& ct1,
    const eConstThermo& ct2
)
{
    specieThermo st(((const specieThermo&)ct1) - ((const specieThermo&)ct2));

    return eConstThermo
    (
        st,
        ct1.nMoles()/st.nMoles()*ct1.Cv_ - ct2.nMoles()/st.nMoles()*ct2.Cv_,
        ct1.nMoles()/st.nMoles()*ct1.Hf_ - ct2.nMoles()/st.nMoles()*ct2.Hf_
    );
}


inline Foam::eConstThermo Foam::operator*
(
    const scalar s,
    const eConstThermo& ct
)
{
    return eConstThermo
    (
        s*((const specieThermo&)ct),
        ct.Cv_,
        ct.Hf_
    );
}


inline Foam::eConstThermo Foam::operator==
(
    const eConstThermo& ct1,
    const eConstThermo& ct2
)
{
    return ct2 - ct1;
}


// ************************************************************************* //
