/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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

#include "polynomialThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline Foam::polynomialThermo<EquationOfState, PolySize>::polynomialThermo
(
    const EquationOfState& pt,
    const scalar Hf,
    const Polynomial<PolySize>& cpPoly,
    const typename Polynomial<PolySize>::intPolyType& dhPoly,
    const Polynomial<PolySize>& sPoly
)
:
    EquationOfState(pt),
    Hf_(Hf),
    cpPolynomial_(cpPoly),
    dhPolynomial_(dhPoly),
    sPolynomial_(sPoly)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline Foam::polynomialThermo<EquationOfState, PolySize>::polynomialThermo
(
    const word& name,
    const polynomialThermo& pt
)
:
    EquationOfState(name, pt),
    Hf_(pt.Hf_),
    cpPolynomial_(pt.cpPolynomial_),
    dhPolynomial_(pt.dhPolynomial_),
    sPolynomial_(pt.sPolynomial_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::polynomialThermo<EquationOfState, PolySize>::cp
(
    const scalar T
) const
{
    return cpPolynomial_.evaluate(T)*this->W();
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::polynomialThermo<EquationOfState, PolySize>::h
(
    const scalar T
) const
{
    return (dhPolynomial_.evaluate(T) + Hf_)*this->W();
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::polynomialThermo<EquationOfState, PolySize>::hs
(
    const scalar T
) const
{
    return dhPolynomial_.evaluate(T)*this->W();
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::polynomialThermo<EquationOfState, PolySize>::hc()
const
{
    return Hf_*this->W();
}


template<class EquationOfState, int PolySize>
inline Foam::scalar Foam::polynomialThermo<EquationOfState, PolySize>::s
(
    const scalar T
) const
{
    WarningIn("scalar polynomialThermo<EquationOfState>::s(const scalar)")
        << "Standard entropy of formation currently set to zero"
        << endl;

    scalar Sf = 0.0;

    return (sPolynomial_.evaluate(T) + Sf)*this->W();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline void Foam::polynomialThermo<EquationOfState, PolySize>::operator+=
(
    const polynomialThermo<EquationOfState, PolySize>& pt
)
{
    scalar molr1 = this->nMoles();

    EquationOfState::operator+=(pt);

    molr1 /= this->nMoles();
    scalar molr2 = pt.nMoles()/this->nMoles();

    Hf_ = molr1*Hf_ + molr2*pt.Hf_;
    cpPolynomial_ = molr1*cpPolynomial_ + molr2*pt.cpPolynomial_;
    dhPolynomial_ = molr1*dhPolynomial_ + molr2*pt.dhPolynomial_;
    sPolynomial_ = molr1*sPolynomial_ + molr2*pt.sPolynomial_;
}


template<class EquationOfState, int PolySize>
inline void Foam::polynomialThermo<EquationOfState, PolySize>::operator-=
(
    const polynomialThermo<EquationOfState, PolySize>& pt
)
{
    scalar molr1 = this->nMoles();

    EquationOfState::operator-=(pt);

    molr1 /= this->nMoles();
    scalar molr2 = pt.nMoles()/this->nMoles();

    Hf_ = molr1*Hf_ - molr2*pt.Hf_;
    cpPolynomial_ = molr1*cpPolynomial_ - molr2*pt.cpPolynomial_;
    dhPolynomial_ = molr1*dhPolynomial_ - molr2*pt.dhPolynomial_;
    sPolynomial_ = molr1*sPolynomial_ - molr2*pt.sPolynomial_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState, int PolySize>
inline Foam::polynomialThermo<EquationOfState, PolySize> Foam::operator+
(
    const polynomialThermo<EquationOfState, PolySize>& pt1,
    const polynomialThermo<EquationOfState, PolySize>& pt2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(pt1)
      + static_cast<const EquationOfState&>(pt2)
    );

    scalar molr1 = pt1.nMoles()/eofs.nMoles();
    scalar molr2 = pt2.nMoles()/eofs.nMoles();
    return polynomialThermo<EquationOfState, PolySize>
    (
        eofs,
        molr1*pt1.Hf_ + molr2*pt2.Hf_,
        molr1*pt1.cpPolynomial_ + molr2*pt2.cpPolynomial_,
        molr1*pt1.dhPolynomial_ + molr2*pt2.dhPolynomial_,
        molr1*pt1.sPolynomial_ + molr2*pt2.sPolynomial_
    );
}


template<class EquationOfState, int PolySize>
inline Foam::polynomialThermo<EquationOfState, PolySize> Foam::operator-
(
    const polynomialThermo<EquationOfState, PolySize>& pt1,
    const polynomialThermo<EquationOfState, PolySize>& pt2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(pt1)
      - static_cast<const EquationOfState&>(pt2)
    );

    scalar molr1 = pt1.nMoles()/eofs.nMoles();
    scalar molr2 = pt2.nMoles()/eofs.nMoles();
    return polynomialThermo<EquationOfState, PolySize>
    (
        eofs,
        molr1*pt1.Hf_ - molr2*pt2.Hf_,
        molr1*pt1.cpPolynomial_ - molr2*pt2.cpPolynomial_,
        molr1*pt1.dhPolynomial_ - molr2*pt2.dhPolynomial_,
        molr1*pt1.sPolynomial_ - molr2*pt2.sPolynomial_
    );
}


template<class EquationOfState, int PolySize>
inline Foam::polynomialThermo<EquationOfState, PolySize> Foam::operator*
(
    const scalar s,
    const polynomialThermo<EquationOfState, PolySize>& pt
)
{
    return polynomialThermo<EquationOfState, PolySize>
    (
        s*static_cast<const EquationOfState&>(pt),
        pt.Hf_,
        pt.cpPolynomial_,
        pt.dhPolynomial_,
        pt.sPolynomial_
    );
}


template<class EquationOfState, int PolySize>
inline Foam::polynomialThermo<EquationOfState, PolySize> Foam::operator==
(
    const polynomialThermo<EquationOfState, PolySize>& pt1,
    const polynomialThermo<EquationOfState, PolySize>& pt2
)
{
    return pt2 - pt1;
}


// ************************************************************************* //
