/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "complex.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::pTraits<Foam::complex>::typeName = "complex";
const char* const Foam::pTraits<Foam::complex>::componentNames[] = {"re", "im"};

const Foam::complex Foam::pTraits<Foam::complex>::zero(0, 0);
const Foam::complex Foam::pTraits<Foam::complex>::one(1, 0);

const Foam::complex Foam::pTraits<Foam::complex>::min(-VGREAT, -VGREAT);
const Foam::complex Foam::pTraits<Foam::complex>::max(VGREAT, VGREAT);

const Foam::complex Foam::pTraits<Foam::complex>::rootMin
(
    -ROOTVGREAT, -ROOTVGREAT
);

const Foam::complex Foam::pTraits<Foam::complex>::rootMax
(
    ROOTVGREAT, ROOTVGREAT
);


Foam::pTraits<Foam::complex>::pTraits(const complex& val)
:
    p_(val)
{}


Foam::pTraits<Foam::complex>::pTraits(Istream& is)
{
    is >> p_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::complex::complex(Istream& is)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::word Foam::name(const complex& c)
{
    return '(' + std::to_string(c.Re()) + ',' + std::to_string(c.Im()) + ')';
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, complex& c)
{
    scalar r, i;

    is.readBegin("complex");
    is >> r >> i;
    is.readEnd("complex");

    c.real(r);
    c.imag(i);

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const complex& c)
{
    os  << token::BEGIN_LIST
        << c.real() << token::SPACE << c.imag()
        << token::END_LIST;

    return os;
}


// ************************************************************************* //
