/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Description
    Dimension set for the base types.
    This type may be used to implement rigorous dimension checking
    for algebraic manipulation.

\*---------------------------------------------------------------------------*/

#include "dimensionSet.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dimensionSet, 1);

const scalar dimensionSet::smallExponent = SMALL;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dimensionSet::dimensionSet
(
    const scalar mass,
    const scalar length,
    const scalar time,
    const scalar temperature,
    const scalar moles,
    const scalar current,
    const scalar luminousIntensity
)
{
    exponents_[MASS] = mass;
    exponents_[LENGTH] = length;
    exponents_[TIME] = time;
    exponents_[TEMPERATURE] = temperature;
    exponents_[MOLES] = moles;
    exponents_[CURRENT] = current;
    exponents_[LUMINOUS_INTENSITY] = luminousIntensity;
}


dimensionSet::dimensionSet
(
    const scalar mass,
    const scalar length,
    const scalar time,
    const scalar temperature,
    const scalar moles
)
{
    exponents_[MASS] = mass;
    exponents_[LENGTH] = length;
    exponents_[TIME] = time;
    exponents_[TEMPERATURE] = temperature;
    exponents_[MOLES] = moles;
    exponents_[CURRENT] = 0;
    exponents_[LUMINOUS_INTENSITY] = 0;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool dimensionSet::dimensionless() const
{
    bool Dimensionless = true;

    for (int Dimension=0; Dimension<nDimensions; Dimension++)
    {
        Dimensionless = Dimensionless &&
        (
            exponents_[Dimension] < smallExponent
         && exponents_[Dimension] > -smallExponent
        );
    }

    return Dimensionless;
}


void dimensionSet::reset(const dimensionSet& ds)
{
    for (int Dimension=0; Dimension<nDimensions; Dimension++)
    {
        exponents_[Dimension] = ds.exponents_[Dimension];
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

scalar dimensionSet::operator[](const dimensionType type) const
{
    return exponents_[type];
}

scalar& dimensionSet::operator[](const dimensionType type)
{
    return exponents_[type];
}


bool dimensionSet::operator==(const dimensionSet& ds) const
{
    bool equall = true;

    for (int Dimension=0; Dimension<nDimensions; Dimension++)
    {
        equall = equall &&
            (mag(exponents_[Dimension] - ds.exponents_[Dimension])
          < smallExponent);
    }

    return equall;
}

bool dimensionSet::operator!=(const dimensionSet& ds) const
{
    return !(operator==(ds));
}


bool dimensionSet::operator=(const dimensionSet& ds) const
{
    if (dimensionSet::debug && *this != ds)
    {
        FatalErrorIn("dimensionSet::operator=(const dimensionSet& ds) const")
            << "Different dimensions for =" << endl
            << "     dimensions : " << *this << " = " << ds << endl
            << abort(FatalError);
    }

    return true;
}


bool dimensionSet::operator+=(const dimensionSet& ds) const
{
    if (dimensionSet::debug && *this != ds)
    {
        FatalErrorIn("dimensionSet::operator+=(const dimensionSet& ds) const")
            << "Different dimensions for +=" << endl
            << "     dimensions : " << *this << " = " << ds << endl
            << abort(FatalError);
    }

    return true;
}

bool dimensionSet::operator-=(const dimensionSet& ds) const
{
    if (dimensionSet::debug && *this != ds)
    {
        FatalErrorIn("dimensionSet::operator-=(const dimensionSet& ds) const")
            << "Different dimensions for -=" << endl
            << "     dimensions : " << *this << " = " << ds << endl
            << abort(FatalError);
    }

    return true;
}

bool dimensionSet::operator*=(const dimensionSet& ds)
{
    reset((*this)*ds);

    return true;
}

bool dimensionSet::operator/=(const dimensionSet& ds)
{
    reset((*this)/ds);

    return true;
}


// * * * * * * * * * * * * * * * Friend functions * * * * * * * * * * * * * * //

dimensionSet max(const dimensionSet& ds1, const dimensionSet& ds2)
{
    if (dimensionSet::debug && ds1 != ds2)
    {
        FatalErrorIn("max(const dimensionSet& ds1, const dimensionSet& ds2)")
            << "Arguments of max have different dimensions" << endl
            << "     dimensions : " << ds1 << " and " << ds2 << endl
            << abort(FatalError);
    }

    return ds1;
}

dimensionSet min(const dimensionSet& ds1, const dimensionSet& ds2)
{
    if (dimensionSet::debug && ds1 != ds2)
    {
        FatalErrorIn("min(const dimensionSet& ds1, const dimensionSet& ds2)")
            << "Arguments of min have different dimensions" << endl
            << "     dimensions : " << ds1 << " and " << ds2 << endl
            << abort(FatalError);
    }

    return ds1;
}


dimensionSet cmptMultiply(const dimensionSet& ds1, const dimensionSet& ds2)
{
    return ds1*ds2;
}


dimensionSet cmptDivide(const dimensionSet& ds1, const dimensionSet& ds2)
{
    return ds1/ds2;
}


dimensionSet pow(const dimensionSet& ds, const scalar p)
{
    dimensionSet dimPow
    (
        ds[dimensionSet::MASS]*p,
        ds[dimensionSet::LENGTH]*p,
        ds[dimensionSet::TIME]*p,
        ds[dimensionSet::TEMPERATURE]*p,
        ds[dimensionSet::MOLES]*p,
        ds[dimensionSet::CURRENT]*p,
        ds[dimensionSet::LUMINOUS_INTENSITY]*p
    );

    return dimPow;
}

dimensionSet pow(const dimensionSet& ds, const dimensionedScalar& dS)
{
    if (dimensionSet::debug && !dS.dimensions().dimensionless())
    {
        FatalErrorIn("pow(const dimensionSet& ds, const dimensionedScalar& dS)")
            << "Exponent of pow are not dimensionless"
            << abort(FatalError);
    }

    dimensionSet dimPow
    (
        ds[dimensionSet::MASS]*dS.value(),
        ds[dimensionSet::LENGTH]*dS.value(),
        ds[dimensionSet::TIME]*dS.value(),
        ds[dimensionSet::TEMPERATURE]*dS.value(),
        ds[dimensionSet::MOLES]*dS.value(),
        ds[dimensionSet::CURRENT]*dS.value(),
        ds[dimensionSet::LUMINOUS_INTENSITY]*dS.value()
    );

    return dimPow;
}

dimensionSet pow(const dimensionedScalar& dS, const dimensionSet& ds)
{
    if
    (
        dimensionSet::debug
     && !dS.dimensions().dimensionless()
     && !ds.dimensionless())
    {
        FatalErrorIn("pow(const dimensionedScalar& dS, const dimensionSet& ds)")
            << "Argument or exponent of pow not dimensionless" << endl
            << abort(FatalError);
    }

    return ds;
}


dimensionSet sqr(const dimensionSet& ds)
{
    return pow(ds, 2);
}

dimensionSet pow3(const dimensionSet& ds)
{
    return pow(ds, 3);
}

dimensionSet pow4(const dimensionSet& ds)
{
    return pow(ds, 4);
}

dimensionSet pow5(const dimensionSet& ds)
{
    return pow(ds, 5);
}

dimensionSet pow6(const dimensionSet& ds)
{
    return pow(ds, 6);
}

dimensionSet sqrt(const dimensionSet& ds)
{
    return pow(ds, 0.5);
}

dimensionSet magSqr(const dimensionSet& ds)
{
    return pow(ds, 2);
}

dimensionSet mag(const dimensionSet& ds)
{
    return ds;
}

dimensionSet sign(const dimensionSet&)
{
    return dimless;
}

dimensionSet pos(const dimensionSet&)
{
    return dimless;
}

dimensionSet neg(const dimensionSet&)
{
    return dimless;
}

dimensionSet inv(const dimensionSet& ds)
{
    return dimless/ds;
}

dimensionSet trans(const dimensionSet& ds)
{
    if (dimensionSet::debug && !ds.dimensionless())
    {
        FatalErrorIn("trans(const dimensionSet& ds)")
            << "Argument of trancendental function not dimensionless"
            << abort(FatalError);
    }

    return ds;
}

dimensionSet transform(const dimensionSet& ds)
{
    return ds;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

dimensionSet operator-(const dimensionSet& ds)
{
    return ds;
}

dimensionSet operator+(const dimensionSet& ds1, const dimensionSet& ds2)
{
    dimensionSet dimSum(ds1);

    if (dimensionSet::debug && ds1 != ds2)
    {
        FatalErrorIn
            ("operator+(const dimensionSet& ds1, const dimensionSet& ds2)")
            << "LHS and RHS of + have different dimensions" << endl
            << "     dimensions : " << ds1 << " + " << ds2 << endl
            << abort(FatalError);
    }

    return dimSum;
}

dimensionSet operator-(const dimensionSet& ds1, const dimensionSet& ds2)
{
    dimensionSet dimDifference(ds1);

    if (dimensionSet::debug && ds1 != ds2)
    {
        FatalErrorIn
            ("operator-(const dimensionSet& ds1, const dimensionSet& ds2)")
            << "LHS and RHS of - have different dimensions" << endl
            << "     dimensions : " << ds1 << " - " << ds2 << endl
            << abort(FatalError);
    }

    return dimDifference;
}

dimensionSet operator*(const dimensionSet& ds1, const dimensionSet& ds2)
{
    dimensionSet dimProduct(ds1);

    for (int Dimension=0; Dimension<dimensionSet::nDimensions; Dimension++)
    {
        dimProduct.exponents_[Dimension] += ds2.exponents_[Dimension];
    }

    return dimProduct;
}

dimensionSet operator/(const dimensionSet& ds1, const dimensionSet& ds2)
{
    dimensionSet dimQuotient(ds1);

    for (int Dimension=0; Dimension<dimensionSet::nDimensions; Dimension++)
    {
        dimQuotient.exponents_[Dimension] -= ds2.exponents_[Dimension];
    }

    return dimQuotient;
}


dimensionSet operator&(const dimensionSet& ds1, const dimensionSet& ds2)
{
    return ds1*ds2;
}

dimensionSet operator^(const dimensionSet& ds1, const dimensionSet& ds2)
{
    return ds1*ds2;
}

dimensionSet operator&&(const dimensionSet& ds1, const dimensionSet& ds2)
{
    return ds1*ds2;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
