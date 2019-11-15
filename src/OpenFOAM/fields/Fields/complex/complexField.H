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

Typedef
    Foam::complexField

Description
    Specialisation of Field\<T\> for complex.

SourceFiles
    complexField.C

\*---------------------------------------------------------------------------*/

#ifndef complexField_H
#define complexField_H

#include "complex.H"
#include "scalarField.H"

#define TEMPLATE
#include "FieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef Field<complex> complexField;

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Zip together complex field from components
void zip
(
    complexField& result,
    const UList<scalar>& re,
    const UList<scalar>& im
);

//- Unzip complex field into components
void unzip
(
    const UList<complex>& input,
    scalarField& re,
    scalarField& im
);


//- Zip up two lists of values into a list of complex
complexField ComplexField
(
    const UList<scalar>& re,
    const UList<scalar>& im
);

//- Create complex field from a list of real (using imag == 0)
complexField ReComplexField(const UList<scalar>& re);

//- Create complex field from a list of imag (using real == 0)
complexField ImComplexField(const UList<scalar>& im);

//- Extract real component
scalarField Re(const UList<complex>& cf);

//- Extract imag component
scalarField Im(const UList<complex>& cf);

//- Sum real and imag components
scalarField ReImSum(const UList<complex>& cf);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Sum product
template<>
complex sumProd(const UList<complex>& f1, const UList<complex>& f2);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

BINARY_TYPE_OPERATOR(complex, complex, complex, +, add)
BINARY_TYPE_OPERATOR(complex, complex, complex, -, subtract)

BINARY_OPERATOR(complex, complex, complex, *, multiply)
BINARY_OPERATOR(complex, complex, complex, /, divide)


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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
