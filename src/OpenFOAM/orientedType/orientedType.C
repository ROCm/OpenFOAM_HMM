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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::orientedType::orientedOption
>
Foam::orientedType::orientedOptionNames
{
    { orientedOption::ORIENTED, "oriented" },
    { orientedOption::UNORIENTED, "unoriented" },
    { orientedOption::UNKNOWN, "unknown" },
};


bool Foam::orientedType::checkType
(
    const orientedType& ot1,
    const orientedType& ot2
)
{
    if
    (
        (ot1.oriented() == UNKNOWN)
     || (ot2.oriented() == UNKNOWN)
     || (ot1.oriented() == ot2.oriented())
    )
    {
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::orientedType::orientedType()
:
    oriented_(UNKNOWN)
{}


Foam::orientedType::orientedType(const orientedType& ot)
:
    oriented_(ot.oriented_)
{}


Foam::orientedType::orientedType(const bool oriented)
:
    oriented_(oriented ? ORIENTED : UNORIENTED)
{}


Foam::orientedType::orientedType(Istream& is)
:
    oriented_(orientedOptionNames.read(is))
{
    is.check(FUNCTION_NAME);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::orientedType::orientedOption& Foam::orientedType::oriented()
{
    return oriented_;
}


Foam::orientedType::orientedOption Foam::orientedType::oriented() const
{
    return oriented_;
}


void Foam::orientedType::setOriented(const bool oriented)
{
    oriented_ = oriented ? ORIENTED : UNORIENTED;
}


void Foam::orientedType::read(const dictionary& dict)
{
    oriented_ = orientedOptionNames.lookupOrDefault
    (
        "oriented",
        dict,
        orientedOption::UNKNOWN
    );
}


void Foam::orientedType::writeEntry(Ostream& os) const
{
    if (oriented_ == ORIENTED)
    {
        os.writeEntry("oriented", orientedOptionNames[oriented_]);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::orientedType::operator=(const orientedType& ot)
{
    // Oriented state is inherited on assignment
    oriented_ = ot.oriented();
}


void Foam::orientedType::operator+=(const orientedType& ot)
{
    // Set the oriented status if it was unknown
    if (oriented_ == UNKNOWN)
    {
        oriented_ = ot.oriented();
    }

    if (!checkType(*this, ot))
    {
        FatalErrorInFunction
            << "Operator += is undefined for "
            << orientedOptionNames[oriented_] << " and "
            << orientedOptionNames[ot.oriented()] << " types"
            << abort(FatalError);
    }
}


void Foam::orientedType::operator-=(const orientedType& ot)
{
    // Set the oriented status if it was unknown
    if (oriented_ == UNKNOWN)
    {
        oriented_ = ot.oriented();
    }

    if (!checkType(*this, ot))
    {
        FatalErrorInFunction
            << "Operator -= is undefined for "
            << orientedOptionNames[oriented_] << " and "
            << orientedOptionNames[ot.oriented()] << " types"
            << abort(FatalError);
    }
}


void Foam::orientedType::operator*=(const orientedType& ot)
{
    const orientedType& ot1 = *this;
    if (ot1() ^ ot())
    {
        oriented_ = ORIENTED;
    }
    else
    {
        oriented_ = UNORIENTED;
    }
}


void Foam::orientedType::operator/=(const orientedType& ot)
{
    const orientedType& ot1 = *this;
    if (ot1() ^ ot())
    {
        oriented_ = ORIENTED;
    }
    else
    {
        oriented_ = UNORIENTED;
    }
}


void Foam::orientedType::operator*=(const scalar s)
{
    // No change to oriented_ flag
}


void Foam::orientedType::operator/=(const scalar s)
{
    // No change to oriented_ flag
}


bool Foam::orientedType::operator()() const
{
    return oriented_ == ORIENTED;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

Foam::orientedType Foam::max(const orientedType& ot1, const orientedType& ot2)
{
    if (!orientedType::checkType(ot1, ot2))
    {
        FatalErrorInFunction
            << "Operator max is undefined for "
            << orientedType::orientedOptionNames[ot1.oriented()] << " and "
            << orientedType::orientedOptionNames[ot2.oriented()] << " types"
            << abort(FatalError);
    }

    return ot1;
}


Foam::orientedType Foam::min(const orientedType& ot1, const orientedType& ot2)
{
    if (!orientedType::checkType(ot1, ot2))
    {
        FatalErrorInFunction
            << "Operator min is undefined for "
            << orientedType::orientedOptionNames[ot1.oriented()] << " and "
            << orientedType::orientedOptionNames[ot2.oriented()] << "types"
            << abort(FatalError);
    }

    return ot1;
}


Foam::orientedType Foam::cmptMultiply
(
    const orientedType& ot1,
    const orientedType& ot2
)
{
    return ot1 ^ ot2;
}


Foam::orientedType Foam::cmptDivide
(
    const orientedType& ot1,
    const orientedType& ot2
)
{
    return ot1 ^ ot2;
}


Foam::orientedType Foam::cmptAv(const orientedType& ot)
{
    return ot;
}


Foam::orientedType Foam::pow(const orientedType& ot, const scalar r)
{
    // Undefined???
    // - only defined for integers where:
    //   - odd powers = oriented_ = yes (if ot is oriented)
    //   - even powers = oriented_ = no
    return ot;
}


Foam::orientedType Foam::sqr(const orientedType& ot)
{
    return orientedType(false);
}


Foam::orientedType Foam::pow3(const orientedType& ot)
{
    return ot;
}


Foam::orientedType Foam::pow4(const orientedType& ot)
{
    return orientedType(false);
}


Foam::orientedType Foam::pow5(const orientedType& ot)
{
    return ot;
}


Foam::orientedType Foam::pow6(const orientedType& ot)
{
    return orientedType(false);
}


Foam::orientedType Foam::pow025(const orientedType& ot)
{
    return ot;
}


Foam::orientedType Foam::sqrt(const orientedType& ot)
{
    return ot;
}


Foam::orientedType Foam::cbrt(const orientedType& ot)
{
    return ot;
}


Foam::orientedType Foam::magSqr(const orientedType& ot)
{
    return orientedType(false);
}


Foam::orientedType  Foam::mag(const orientedType& ot)
{
    return orientedType(false);
}


Foam::orientedType  Foam::sign(const orientedType& ot)
{
    return ot;
}


Foam::orientedType  Foam::pos(const orientedType& ot)
{
    return ot;
}


Foam::orientedType  Foam::pos0(const orientedType& ot)
{
    return ot;
}


Foam::orientedType  Foam::neg(const orientedType& ot)
{
    return ot;
}


Foam::orientedType  Foam::neg0(const orientedType& ot)
{
    return ot;
}


Foam::orientedType  Foam::posPart(const orientedType& ot)
{
    return ot;
}


Foam::orientedType  Foam::negPart(const orientedType& ot)
{
    return ot;
}


Foam::orientedType  Foam::inv(const orientedType& ot)
{
    return ot;
}


Foam::orientedType Foam::trans(const orientedType& ot)
{
    return ot;
}


Foam::orientedType Foam::atan2
(
    const orientedType& ot1,
    const orientedType& ot2
)
{
    if (!orientedType::checkType(ot1, ot2))
    {
        FatalErrorInFunction
            << "Operator atan2 is undefined for "
            << orientedType::orientedOptionNames[ot1.oriented()] << " and "
            << orientedType::orientedOptionNames[ot2.oriented()] << "types"
            << abort(FatalError);
    }

    return ot1;
}


Foam::orientedType Foam::transform(const orientedType& ot)
{
    return ot;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, orientedType& ot)
{
    ot.oriented_ = orientedType::orientedOptionNames.read(is);

    is.check(FUNCTION_NAME);

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const orientedType& ot)
{
    os << orientedType::orientedOptionNames[ot.oriented()];

    os.check(FUNCTION_NAME);

    return os;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

Foam::orientedType Foam::operator+
(
    const orientedType& ot1,
    const orientedType& ot2
)
{
    if (!orientedType::checkType(ot1, ot2))
    {
        FatalErrorInFunction
            << "Operator + is undefined for "
            << orientedType::orientedOptionNames[ot1.oriented()] << " and "
            << orientedType::orientedOptionNames[ot2.oriented()] << " types"
            << abort(FatalError);
    }

    // Note use of () operators to convert to boolean op
    return orientedType(ot1() || ot2());
}


Foam::orientedType Foam::operator-(const orientedType& ot)
{
    return ot;
}


Foam::orientedType Foam::operator-
(
    const orientedType& ot1,
    const orientedType& ot2
)
{
    if (!orientedType::checkType(ot1, ot2))
    {
        FatalErrorInFunction
            << "Operator - is undefined for "
            << orientedType::orientedOptionNames[ot1.oriented()] << " and "
            << orientedType::orientedOptionNames[ot2.oriented()] << " types"
            << abort(FatalError);
    }

    // Note use of () operators to convert to boolean op
    return orientedType(ot1() || ot2());
}


Foam::orientedType Foam::operator*(const scalar s, const orientedType& ot)
{
    return ot;
}


Foam::orientedType Foam::operator/(const orientedType& ot, const scalar s)
{
    return ot;
}


Foam::orientedType Foam::operator/
(
    const orientedType& ot1,
    const orientedType& ot2
)
{
    return ot1 ^ ot2;
}


Foam::orientedType Foam::operator*
(
    const orientedType& ot1,
    const orientedType& ot2
)
{
    return ot1 ^ ot2;
}


Foam::orientedType Foam::operator^
(
    const orientedType& ot1,
    const orientedType& ot2
)
{
    // Note use of () operators to convert to boolean op
    return orientedType(ot1() ^ ot2());
}


Foam::orientedType Foam::operator&
(
    const orientedType& ot1,
    const orientedType& ot2
)
{
    return ot1 ^ ot2;
}


Foam::orientedType Foam::operator&&
(
    const orientedType& ot1,
    const orientedType& ot2
)
{
    return orientedType(false);
}


// ************************************************************************* //
