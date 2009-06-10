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

#include "Polynomial.H"
#include "error.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<int PolySize>
Foam::Polynomial<PolySize>::Polynomial()
:
    VectorSpace<Polynomial<PolySize>, scalar, PolySize>(),
    name_("unknownPolynomialName")
{}


template<int PolySize>
Foam::Polynomial<PolySize>::Polynomial(const word& name, Istream& is)
:
    VectorSpace<Polynomial<PolySize>, scalar, PolySize>(),
    name_(is)
{
    if (name_ != name)
    {
        FatalErrorIn
        (
            "Foam::Polynomial<PolySize>::Polynomial(const word&, Istream&)"
        )   << "Expected polynomial name " << name << " but read " << name_
            << nl << exit(FatalError);
    }

    VectorSpace<Polynomial<PolySize>, scalar, PolySize>::
        operator=(polyType(is));

    if (this->size() == 0)
    {
        FatalErrorIn
        (
            "Foam::Polynomial<PolySize>::Polynomial(const word&, Istream&)"
        )   << "Polynomial coefficients for entry " << name_
            << " are invalid (empty)" << nl << exit(FatalError);
    }
}


template<int PolySize>
Foam::Polynomial<PolySize>::Polynomial(const Polynomial<PolySize>& poly)
:
    VectorSpace<Polynomial<PolySize>, scalar, PolySize>(poly),
    name_(poly.name_)
{}


template<int PolySize>
Foam::Polynomial<PolySize>::Polynomial
(
    const word& name,
    const Polynomial<PolySize>& poly
)
:
    VectorSpace<Polynomial<PolySize>, scalar, PolySize>(poly),
    name_(name)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<int PolySize>
Foam::Polynomial<PolySize>::~Polynomial()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<int PolySize>
const Foam::word& Foam::Polynomial<PolySize>::name() const
{
    return name_;
}


template<int PolySize>
Foam::scalar Foam::Polynomial<PolySize>::evaluate(const scalar x) const
{
    scalar y = 0.0;
    forAll(*this, i)
    {
        y += this->v_[i]*pow(x, i);
    }

    return y;
}


template<int PolySize>
Foam::scalar Foam::Polynomial<PolySize>::integrateLimits
(
    const scalar x1,
    const scalar x2
) const
{
    scalar intx = 0.0;

    forAll(*this, i)
    {
        intx += this->v_[i]/(i + 1)*(pow(x2, i + 1) - pow(x1, i + 1));
    }

    return intx;
}


template<int PolySize>
typename Foam::Polynomial<PolySize>::intPolyType
Foam::Polynomial<PolySize>::integrate(const scalar intConstant)
{
    intPolyType newCoeffs;

    newCoeffs[0] = intConstant;
    forAll(*this, i)
    {
        newCoeffs[i + 1] = this->v_[i]/(i + 1);
    }

    return newCoeffs;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<int PolySize>
void Foam::Polynomial<PolySize>::operator=(const Polynomial<PolySize>& poly)
{
    name_ = poly.name_;
    VectorSpace<Polynomial<PolySize>, scalar, PolySize>::operator=(poly);
}


// ************************************************************************* //
