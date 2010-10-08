/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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
    PolynomialTest

Description
    Test application for the templated Polynomial class

\*---------------------------------------------------------------------------*/

#include "IFstream.H"
#include "Polynomial.H"
#include "Random.H"
#include "cpuTime.H"

using namespace Foam;

const int nCoeffs = 8;
const scalar coeff[] = { 0.11, 0.45, -0.94, 1.58, 2.58, 0.08, 3.15, -4.78 };


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


scalar polyValue1(const scalar x)
{
    // "normal" evaluation using pow()
    scalar value = coeff[0];

    for (int i=1; i < nCoeffs; ++i)
    {
        value += coeff[i]*pow(x, i);
    }

    return value;
}


// calculation avoiding pow()
scalar polyValue2(const scalar x)
{
    scalar value = coeff[0];

    scalar powX = x;
    for (int i=1; i < nCoeffs; ++i, powX *= x)
    {
        value += coeff[i] * powX;
    }

    return value;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    const label n = 10000;
    const label nIters = 1000;
    scalar sum = 0.0;

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


    //
    // test speed of Polynomial:
    //
    Info<< "start timing loops" << nl
        << "~~~~~~~~~~~~~~~~~~" << endl;

    cpuTime timer;

    for (int loop = 0; loop < n; ++loop)
    {
        sum = 0.0;
        for (label iter = 0; iter < nIters; ++iter)
        {
            sum += poly.evaluate(loop+iter);
        }
    }
    Info<< "evaluate:     " << sum
        << " in " << timer.cpuTimeIncrement() << " s\n";


    for (int loop = 0; loop < n; ++loop)
    {
        sum = 0.0;
        for (label iter = 0; iter < nIters; ++iter)
        {
            sum += polyValue(loop+iter);
        }
    }
    Info<< "hard-coded 0: " << sum
        << " in " << timer.cpuTimeIncrement() << " s\n";


    for (int loop = 0; loop < n; ++loop)
    {
        sum = 0.0;
        for (label iter = 0; iter < nIters; ++iter)
        {
            sum += polyValue1(loop+iter);
        }
    }
    Info<< "hard-coded 1: " << sum
        << " in " << timer.cpuTimeIncrement() << " s\n";

    for (int loop = 0; loop < n; ++loop)
    {
        sum = 0.0;
        for (label iter = 0; iter < nIters; ++iter)
        {
            sum += polyValue2(loop+iter);
        }
    }
    Info<< "hard-coded 2: " << sum
        << " in " << timer.cpuTimeIncrement() << " s\n";


    Info<< nl << "Done." << endl;

    return 0;
}


// ************************************************************************* //
