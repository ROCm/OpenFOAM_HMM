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
    Some tests for complex numbers
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "complex.H"
#include "complexVector.H"
#include "Field.H"

using namespace Foam;

void print1(const complex& z)
{
    Info<<"r: " << z.real() << " i: " << z.imag() << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    Info<< "complex()     : " << complex() << nl
        << "complex(zero) : " << complex(Zero) << nl
        << "complex::zero : " << complex::zero << nl
        << "complex::one  : " << complex::one << nl
        << "complex(scalar) : " << complex(3.14519) << nl
        << nl;

    Info<< "complexVector::zero : " << complexVector::zero << nl
        << "complexVector::one  : " << complexVector::one << nl
        << nl;

    for (complex c : { complex{1, 0}, complex{1, 2}} )
    {
        print1(c);

        // TDB: allow implicit construct from scalar?
        //
        // if (c == 1.0)
        // {
        //     Info<< c << " == " << 1 << nl;
        // }
    }

    Field<complex> fld1(3, complex(2.0, 1.0));

    Info<< "Field " << flatOutput(fld1) << nl;

    Info<< "sum = " << sum(fld1) << nl;
    // Not yet Info<< "min = " << min(fld1) << nl;

    fld1 *= 10;
    Info<< "Multiply: " << flatOutput(fld1) << nl;

    for (complex& c : fld1)
    {
        c = ~c;
    }

    Info<< "Conjugate: " << flatOutput(fld1) << nl;

    Info<< "\nEnd\n" << endl;
    return 0;
}

// ************************************************************************* //
