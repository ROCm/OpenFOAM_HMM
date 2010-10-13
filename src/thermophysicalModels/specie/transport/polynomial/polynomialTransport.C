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

#include "polynomialTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
Foam::polynomialTransport<Thermo, PolySize>::polynomialTransport(Istream& is)
:
    Thermo(is),
    muPolynomial_("muPolynomial", is),
    kappaPolynomial_("kappaPolynomial", is)
{
    muPolynomial_ *= this->W();
    kappaPolynomial_ *= this->W();
}


template<class Thermo, int PolySize>
Foam::polynomialTransport<Thermo, PolySize>::polynomialTransport
(
    const dictionary& dict
)
:
    Thermo(dict),
    muPolynomial_(dict.lookup("muPolynomial")),
    kappaPolynomial_(dict.lookup("kappaPolynomial"))
{
    muPolynomial_ *= this->W();
    kappaPolynomial_ *= this->W();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
void Foam::polynomialTransport<Thermo, PolySize>::write(Ostream& os) const
{
    os  << this->name() << endl;
    os  << token::BEGIN_BLOCK << incrIndent << nl;
    Thermo::write(os);
    os.writeKeyword("muPolynomial") << muPolynomial_/this->W()
        << token::END_STATEMENT << nl;
    os.writeKeyword("kappaPolynomial") << kappaPolynomial_/this->W()
        << token::END_STATEMENT << nl;
    os   << decrIndent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo, int PolySize>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const polynomialTransport<Thermo, PolySize>& pt
)
{
    os  << static_cast<const Thermo&>(pt) << tab
        << "muPolynomial" << tab << pt.muPolynomial_/pt.W() << tab
        << "kappaPolynomial" << tab << pt.kappaPolynomial_/pt.W();

    os.check
    (
        "Ostream& operator<<"
        "("
            "Ostream&, "
            "const polynomialTransport<Thermo, PolySize>&"
        ")"
    );

    return os;
}


// ************************************************************************* //
