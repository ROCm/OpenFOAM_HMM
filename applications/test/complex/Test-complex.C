/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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

Application

Description
    Tests for complex numbers

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "complex.H"
#include "complexFields.H"
#include "ops.H"
#include "ListOps.H"

using namespace Foam;

void print1(const complex& z)
{
    Info<<"r: " << z.real() << " i: " << z.imag() << nl;
}


// * * * * * * * * * * * * * * * Main Program  * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Info<< "complex()     : " << complex() << nl
        << "complex(zero) : " << complex(Zero) << nl
        << "pTraits<complex>::zero : " << pTraits<complex>::zero << nl
        << "pTraits<complex>::one  : " << pTraits<complex>::one << nl
        << "complex(scalar) : " << complex(3.14519) << nl
        << nl;

    std::complex<scalar> c1(10, -3);
    Info<< "std::complex : " << c1 << nl;
    Info<< "sin: " << std::sin(c1) << nl;

    Info<< "complexVector::zero : " << complexVector::zero << nl
        << "complexVector::one  : " << complexVector::one << nl
        << nl;

    for (complex c : { complex{1, 0}, complex{1, 2}} )
    {
        Info<< nl;
        print1(c);

        Info<< "sin: " << sin(c) << nl;
        Info<< "pow(3): " << pow(c, 3) << nl;
        Info<< "pow3: " << pow3(c) << nl;
        Info<< "log: " << log(c) << nl;
        Info<< "pow025: " << pow025(c) << nl;

        // TDB: allow implicit construct from scalar?
        //
        // if (c == 1.0)
        // {
        //     Info<< c << " == " << 1 << nl;
        // }
    }

    complexField fld1(3, complex(2.0, 1.0));
    complexField fld2(fld1);

    for (complex& c : fld2)
    {
        c = ~c;
    }

    Info<< "Field " << flatOutput(fld1) << nl;
    Info<< "Conjugate: " << flatOutput(fld2) << nl;

    // Some arbitrary change
    for (complex& c : fld2)
    {
        c.Im() *= 5;
    }


    fld1 *= 10;
    Info<< "scalar multiply: " << flatOutput(fld1) << nl;

    fld1 /= 10;
    Info<< "scalar divide: " << flatOutput(fld1) << nl;

    Info<< "sin: " << sin(fld1) << nl;

    Info<< "operator + : " << (fld1 + fld2) << nl;

    // Some operators are still incomplete

    // Info<< "operator * : " << (fld1 * fld2) << nl;
    // Info<< "operator / : " << (fld1 / fld2) << nl;
    // Info<< "operator / : " << (fld1 / 2) << nl;
    // Info<< "operator / : " << (fld1 / fld2) << nl;
    // Info<< "sqrt   : " << sqrt(fld1) << nl;
    // Info<< "pow(2) : " << pow(fld1, 2) << nl;

    Info<< nl << "Elementary complex arithmetic operations:" << nl;
    {
        complex a(6, 1);
        complex b = a;

        Info << "Compound assignment operations:" << nl;

        // Multiplication
        b *= a;
        Info<< "b *= a:" << tab << "b=" << b << nl;

        // Addition
        b += a;
        Info<< "b += a:" << tab << "b=" << b << nl;

        // Subtraction
        b -= a;
        Info<< "b -= a:" << tab << "b=" << b << nl;

        // Division
        b /= a;
        Info<< "b /= a:" << tab << "b=" << b << nl;

        Info << "Operations with scalars:" << nl;

        Info<< "b=" << b << nl;

        // Scalar multiplication
        b *= 2.0;
        Info<< "b*2 (elementwise multiplication):" << tab << b << nl;

        // Scalar addition
        b += 1.0;
        Info<< "b + 1 (only real part):" << tab << b << nl;

        // Scalar subtraction
        b -= 1.0;
        Info<< "b - 1 (only real part):" << tab << b << nl;

        // Scalar division
        b = 1.0/b;
        Info<< "1/b (elementwise division):" << tab << b << nl;
    }

    Info<< nl << "Other mathematical expressions:" << nl;
    {
        complex a(4.3, -3.14);
        complex b(0, -4.3);

        Info<< "a=" << a << tab << "b=" << b << nl;

        // Square-root
        //Info<< "sqrt(a)=" << sqrt(a) << tab << "sqrt(b)=" << sqrt(b) << nl;

        // Square
        Info<< "sqr(a)=" << sqr(a) << tab << "sqr(b)=" << sqr(b) << nl;

        // n^th power
        //Info<< "pow(a,-1)=" << pow(a,-1) << tab
        //      << "pow(b,-1)=" << pow(b,-1) << nl;

        // Exponential
        //Info<< "exp(a)=" << exp(a) << tab << "exp(b)=" << exp(b) << nl;

        // Natural logarithm
        //Info<< "log(a)=" << log(a) << tab << "log(b)=" << log(b) << nl;
    }

    Info<< nl << "End" << nl;

    // Make some changes
    {
        label i = 1;
        for (complex& c : fld1)
        {
            c.Re() += i;
            c.Im() -= 10 - i;
            ++i;
        }
    }

    Info<< nl
        << "field = " << fld1 << nl;

    Info<< "magSqr = "
        << ListOps::create<scalar>
           (
               fld1,
               [](const complex& c) { return magSqr(c); }
           )
        << nl;

    Info
        << "sum = " << sum(fld1) << nl
        << "min = " << min(fld1) << nl
        << "max = " << max(fld1) << nl;

    // MinMax fails since there is no less comparison operator
    // Info<< "min/max = " << MinMax<complex>(fld1) << nl;

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
