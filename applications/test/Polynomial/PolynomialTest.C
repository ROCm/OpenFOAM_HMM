/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

Application
    PolynomialTest

Description
    Test application for the templated Polynomial class

\*---------------------------------------------------------------------------*/

#include "IFstream.H"
#include "Polynomial.H"
#include "Random.H"

using namespace Foam;

scalar polyValue(const scalar x)
{
    // Hard-coded polynomial 8 coeff (7th order)
    return
        0.11
      + 0.45*x
      - 0.94*sqr(x)
      + 1.58*pow3(x)
      - 2.58*pow4(x)
      + 0.08*pow5(x)
      + 3.15*pow6(x)
      - 4.78*x*pow6(x);
}


scalar intPolyValue(const scalar x)
{
    // Hard-coded integrated form of above polynomial
    return
        0.11*x
      + 0.45/2.0*sqr(x)
      - 0.94/3.0*pow3(x)
      + 1.58/4.0*pow4(x)
      - 2.58/5.0*pow5(x)
      + 0.08/6.0*pow6(x)
      + 3.15/7.0*x*pow6(x)
      - 4.78/8.0*x*x*pow6(x);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    IFstream is("polyTestInput");

    Polynomial<8> poly("testPoly", is);
    Polynomial<9> intPoly(poly.integrate(0.0));

    Info<< "poly = " << poly << endl;
    Info<< "intPoly = " << intPoly << nl << endl;

    Info<< "2*poly = " << 2*poly << endl;
    Info<< "poly+poly = " << poly + poly << nl << endl;

    Info<< "3*poly = " << 3*poly << endl;
    Info<< "poly+poly+poly = " << poly + poly + poly << nl << endl;

    Info<< "3*poly - 2*poly = " << 3*poly - 2*poly << nl << endl;

    Polynomial<8> polyCopy = poly;
    Info<< "poly, polyCopy = " << poly << ", " << polyCopy << nl << endl;
    polyCopy = 2.5*poly;
    Info<< "2.5*polyCopy = " << polyCopy << nl << endl;

    Random rnd(123456);
    for (int i=0; i<10; i++)
    {
        scalar x = rnd.scalar01()*100;

        scalar px = polyValue(x);
        scalar ipx = intPolyValue(x);

        scalar pxTest = poly.evaluate(x);
        scalar ipxTest = intPoly.evaluate(x);

        Info<<"\nx = " << x << endl;
        Info<< "    px, pxTest = " << px << ", " << pxTest << endl;
        Info<< "    ipx, ipxTest = " << ipx << ", " << ipxTest << endl;

        if (mag(px - pxTest) > SMALL)
        {
            Info<< "    *** WARNING: px != pxTest: " << px - pxTest << endl;
        }

        if (mag(ipx - ipxTest) > SMALL)
        {
            Info<< "    *** WARNING: ipx != ipxTest: " << ipx - ipxTest << endl;
        }

        Info<< endl;
    }

    Info<< nl << "Done." << endl;

    return 0;
}


// ************************************************************************* //
