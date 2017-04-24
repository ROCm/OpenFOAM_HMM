/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

#include "orientedType.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::orientedType::orientedType()
:
    oriented_(false)
{}


Foam::orientedType::orientedType(const orientedType& of)
:
    oriented_(of.oriented_)
{}


Foam::orientedType::orientedType(const bool oriented)
:
    oriented_(oriented)
{}


Foam::orientedType::orientedType(Istream& is)
:
    oriented_(readBool(is))
{
    is.check(FUNCTION_NAME);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool& Foam::orientedType::oriented()
{
    return oriented_;
}


bool Foam::orientedType::oriented() const
{
    return oriented_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::orientedType::operator=(const orientedType& of)
{
    oriented_ = of.oriented();
}


void Foam::orientedType::operator+=(const orientedType& of)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    if (oriented_ != of.oriented())
    {
        FatalErrorInFunction
            << "Operator += is undefined for oriented and unoriented types. "
            << "oriented:" << oriented_ << ", of:" << of.oriented()
            << abort(FatalError);
    }

    // No change to oriented_ flag
}


void Foam::orientedType::operator-=(const orientedType& of)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    if (oriented_ != of.oriented())
    {
        FatalErrorInFunction
            << "Operator -= is undefined for oriented and unoriented types. "
            << "oriented:" << oriented_ << ", of:" << of.oriented()
            << abort(FatalError);
    }

    // No change to oriented_ flag
}


void Foam::orientedType::operator*=(const orientedType& of)
{
    oriented_ = oriented_ ^ of.oriented();
}


void Foam::orientedType::operator/=(const orientedType& of)
{
    oriented_ = oriented_ ^ of.oriented();
}


void Foam::orientedType::operator*=(const scalar s)
{
//InfoInFunction << "oriented_: " << oriented_ << endl;
    // No change to oriented_ flag
}


void Foam::orientedType::operator/=(const scalar s)
{
//InfoInFunction << "oriented_: " << oriented_ << endl;
    // No change to oriented_ flag
}


bool Foam::orientedType::operator()() const
{
//InfoInFunction << "oriented_: " << oriented_ << endl;
    return oriented_;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

Foam::orientedType Foam::max(const orientedType& of1, const orientedType& of2)
{
//InfoInFunction << "of1:" << of1.oriented() << ", of2:" << of2.oriented() << endl;
    if (of1.oriented() != of2.oriented())
    {
        FatalErrorInFunction
            << "max is undefined for oriented and unoriented types. "
            << "of1:" << of1.oriented() << ", of2:" << of2.oriented()
            << abort(FatalError);
    }

    return of1;
}


Foam::orientedType Foam::min(const orientedType& of1, const orientedType& of2)
{
//InfoInFunction << "of1:" << of1.oriented() << ", of2:" << of2.oriented() << endl;
    if (of1.oriented() != of2.oriented())
    {
        FatalErrorInFunction
            << "min is undefined for oriented and unoriented types. "
            << "of1:" << of1.oriented() << ", of2:" << of2.oriented()
            << abort(FatalError);
    }

    return of1;
}


Foam::orientedType Foam::cmptMultiply
(
    const orientedType& of1,
    const orientedType& of2
)
{
//InfoInFunction << "of1:" << of1.oriented() << ", of2:" << of2.oriented() << endl;
    return orientedType(of1.oriented() ^ of2.oriented());
}


Foam::orientedType Foam::cmptDivide
(
    const orientedType& of1,
    const orientedType& of2
)
{
//InfoInFunction << "of1:" << of1.oriented() << ", of2:" << of2.oriented() << endl;
    return orientedType(of1.oriented() ^ of2.oriented());
}


Foam::orientedType Foam::cmptAv(const orientedType& of)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    return of;
}


Foam::orientedType Foam::pow(const orientedType& of, const scalar r)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    // Undefined???
    // - only defined for integers where:
    //   - odd powers = oriented_ = yes (if of is oriented)
    //   - even powers = oriented_ = no
    return of;
}


Foam::orientedType Foam::sqr(const orientedType& of)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    return orientedType(false);
}


Foam::orientedType Foam::pow3(const orientedType& of)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    return orientedType(of.oriented());
}


Foam::orientedType Foam::pow4(const orientedType& of)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    return orientedType(false);
}


Foam::orientedType Foam::pow5(const orientedType& of)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    return orientedType(of.oriented());
}


Foam::orientedType Foam::pow6(const orientedType& of)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    return orientedType(false);
}


Foam::orientedType Foam::pow025(const orientedType& of)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    return orientedType(of.oriented());
}


Foam::orientedType Foam::sqrt(const orientedType& of)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    return of;
}


Foam::orientedType Foam::cbrt(const orientedType& of)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    return of;
}


Foam::orientedType Foam::magSqr(const orientedType& of)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    return orientedType(false);
}


Foam::orientedType  Foam::mag(const orientedType& of)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    return orientedType(false);
}


Foam::orientedType  Foam::sign(const orientedType& of)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    return of;
}


Foam::orientedType  Foam::pos(const orientedType& of)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    return of;
}


Foam::orientedType  Foam::neg(const orientedType& of)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    return of;
}


Foam::orientedType  Foam::posPart(const orientedType& of)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    return of;
}


Foam::orientedType  Foam::negPart(const orientedType& of)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    return of;
}


Foam::orientedType  Foam::inv(const orientedType& of)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    return of;
}


Foam::orientedType Foam::trans(const orientedType& of)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    return of;
}


Foam::orientedType Foam::atan2
(
    const orientedType& of1,
    const orientedType& of2
)
{
//InfoInFunction << "of1:" << of1.oriented() << ", of2:" << of2.oriented() << endl;
    if (of1.oriented() != of2.oriented())
    {
        FatalErrorInFunction
            << "atan2 is undefined for oriented and unoriented types. "
            << "of1:" << of1.oriented() << ", of2:" << of2.oriented()
            << abort(FatalError);
    }

    return of1;
}


Foam::orientedType Foam::transform(const orientedType& of)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    return of;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, orientedType& of)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    is >> of.oriented_;

    is.check(FUNCTION_NAME);

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const orientedType& of)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    os << of.oriented();

    os.check(FUNCTION_NAME);

    return os;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

Foam::orientedType Foam::operator+
(
    const orientedType& of1,
    const orientedType& of2
)
{
//InfoInFunction << "of1:" << of1.oriented() << ", of2:" << of2.oriented() << endl;
    if (of1.oriented() != of2.oriented())
    {
        FatalErrorInFunction
            << "Operator + is undefined for oriented and unoriented types. "
            << "of1:" << of1.oriented() << ", of2:" << of2.oriented()
            << abort(FatalError);
    }

    return orientedType(of1.oriented() || of2.oriented());
}


Foam::orientedType Foam::operator-(const orientedType& of)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    return orientedType(of);
}


Foam::orientedType Foam::operator-
(
    const orientedType& of1,
    const orientedType& of2
)
{
//InfoInFunction << "of1:" << of1.oriented() << ", of2:" << of2.oriented() << endl;
    if (of1.oriented() != of2.oriented())
    {
        FatalErrorInFunction
            << "Operator - is undefined for oriented and unoriented types. "
            << "of1:" << of1.oriented() << ", of2:" << of2.oriented()
            << abort(FatalError);
    }

    return orientedType(of1.oriented() || of2.oriented());
}


Foam::orientedType Foam::operator*(const scalar s, const orientedType& of)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    return orientedType(of);
}


Foam::orientedType Foam::operator/(const orientedType& of, const scalar s)
{
//InfoInFunction << "of:" << of.oriented() << endl;
    return orientedType(of);
}


Foam::orientedType Foam::operator/
(
    const orientedType& of1,
    const orientedType& of2
)
{
//InfoInFunction << "of1:" << of1.oriented() << ", of2:" << of2.oriented() << endl;
    return orientedType(of1.oriented() ^ of2.oriented());
}


Foam::orientedType Foam::operator*
(
    const orientedType& of1,
    const orientedType& of2
)
{
//InfoInFunction << "of1:" << of1.oriented() << ", of2:" << of2.oriented() << endl;
    return orientedType(of1.oriented() ^ of2.oriented());
}


Foam::orientedType Foam::operator^
(
    const orientedType& of1,
    const orientedType& of2
)
{
//InfoInFunction << "of1:" << of1.oriented() << ", of2:" << of2.oriented() << endl;
    return orientedType(of1.oriented() ^ of2.oriented());
}


Foam::orientedType Foam::operator&
(
    const orientedType& of1,
    const orientedType& of2
)
{
//InfoInFunction << "of1:" << of1.oriented() << ", of2:" << of2.oriented() << endl;
    return orientedType(of1.oriented() ^ of2.oriented());
}


Foam::orientedType Foam::operator&&
(
    const orientedType& of1,
    const orientedType& of2
)
{
//InfoInFunction << "of1:" << of1.oriented() << ", of2:" << of2.oriented() << endl;
    return orientedType(false);
}


// ************************************************************************* //
