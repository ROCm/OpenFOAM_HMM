/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010, 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011 OpenFOAM Foundation
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

#include "complexField.H"
#include "addToRunTimeSelectionTable.H"

#define TEMPLATE
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineCompoundTypeName(List<complex>, complexList);
    addCompoundToRunTimeSelectionTable(List<complex>, complexList);
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::complexField Foam::ComplexField
(
    const UList<scalar>& re,
    const UList<scalar>& im
)
{
    complexField cf(re.size());

    forAll(cf, i)
    {
        cf[i].Re() = re[i];
        cf[i].Im() = im[i];
    }

    return cf;
}


Foam::complexField Foam::ReComplexField(const UList<scalar>& re)
{
    complexField cf(re.size());

    forAll(cf, i)
    {
        cf[i].Re() = re[i];
        cf[i].Im() = 0.0;
    }

    return cf;
}


Foam::complexField Foam::ImComplexField(const UList<scalar>& im)
{
    complexField cf(im.size());

    forAll(cf, i)
    {
        cf[i].Re() = 0.0;
        cf[i].Im() = im[i];
    }

    return cf;
}


Foam::scalarField Foam::ReImSum(const UList<complex>& cf)
{
    scalarField sf(cf.size());

    forAll(sf, i)
    {
        sf[i] = cf[i].Re() + cf[i].Im();
    }

    return sf;
}


Foam::scalarField Foam::Re(const UList<complex>& cf)
{
    scalarField sf(cf.size());

    forAll(sf, i)
    {
        sf[i] = cf[i].Re();
    }

    return sf;
}


Foam::scalarField Foam::Im(const UList<complex>& cf)
{
    scalarField sf(cf.size());

    forAll(sf, i)
    {
        sf[i] = cf[i].Im();
    }

    return sf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
complex sumProd(const UList<complex>& f1, const UList<complex>& f2)
{
    complex result = Zero;
    if (f1.size() && (f1.size() == f2.size()))
    {
        TFOR_ALL_S_OP_F_OP_F(complex, result, +=, complex, f1, *, complex, f2)
    }
    return result;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

UNARY_FUNCTION(complex, complex, pow3)
UNARY_FUNCTION(complex, complex, pow4)
UNARY_FUNCTION(complex, complex, pow5)
UNARY_FUNCTION(complex, complex, pow6)
UNARY_FUNCTION(complex, complex, pow025)
UNARY_FUNCTION(complex, complex, sqrt)
UNARY_FUNCTION(complex, complex, exp)
UNARY_FUNCTION(complex, complex, log)
UNARY_FUNCTION(complex, complex, log10)
UNARY_FUNCTION(complex, complex, sin)
UNARY_FUNCTION(complex, complex, cos)
UNARY_FUNCTION(complex, complex, tan)
UNARY_FUNCTION(complex, complex, asin)
UNARY_FUNCTION(complex, complex, acos)
UNARY_FUNCTION(complex, complex, atan)
UNARY_FUNCTION(complex, complex, sinh)
UNARY_FUNCTION(complex, complex, cosh)
UNARY_FUNCTION(complex, complex, tanh)
UNARY_FUNCTION(complex, complex, asinh)
UNARY_FUNCTION(complex, complex, acosh)
UNARY_FUNCTION(complex, complex, atanh)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
