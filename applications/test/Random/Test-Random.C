/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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
    Test-Random

Description
    Simple test for sequence of random numbers

\*---------------------------------------------------------------------------*/

#include "Rand48.H"
#include "Random.H"

#include <cstdlib>
#include <iostream>
#include <iomanip>

#define TEST_RAW_IEEE

// Construct a positive double with the 48 random bits distributed over
// its fractional part so the resulting FP number is [0.0,1.0).
//
// As per glibc erand48() implementation

#ifdef TEST_RAW_IEEE
#include <ieee754.h>
double randomFraction(const uint64_t bits)
{
    // 48-bit value
    unsigned short int xsubi[3];
    xsubi[0] = (bits & 0xffff);
    xsubi[1] = ((bits >> 16) & 0xffff);
    xsubi[2] = ((bits >> 32) & 0xffff);

    union ieee754_double temp;

    temp.ieee.negative = 0;
    temp.ieee.exponent = IEEE754_DOUBLE_BIAS;
    temp.ieee.mantissa0 = (xsubi[2] << 4) | (xsubi[1] >> 12);
    temp.ieee.mantissa1 = ((xsubi[1] & 0xfff) << 20) | (xsubi[0] << 4);

    // The lower 4 bits of mantissa1 are always 0.
    return temp.d - 1.0;
}
#endif


using namespace Foam;

// Test uniformity of random
void testPosition(const label n)
{
    List<label> samples(n, Zero);

    Random rnd(123456);
    for (label i=0; i < 100000*n; ++i)
    {
        ++samples[rnd.position<label>(0,n-1)];
    }

    Info<< nl << "uniform [0," << n << ")\n  "
        << flatOutput(samples) << nl;
}


// Output with cout instead of Info to keep values unsigned on output
using std::cout;
using std::setw;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Allow multiple passes to ensure reset is working properly
    const label maxIter = 1;

    std::default_random_engine deflt(123456);
    std::mt19937 mtwist(123456);
    std::uniform_real_distribution<scalar> uniform01;

    {
        Rand48 rnd(123456);

        for (label iter=0; iter < maxIter; ++iter)
        {
            rnd.seed(123456);
            ::srand48(123456);
            mtwist.seed(123456);

            Info<< nl << "32-bit random with seed = 123456" << nl;

            cout<< setw(12) << "Rand48()"
                << setw(12) << "lrand48()"
                << setw(12) << "mtwister"
                << setw(12) << "default" << nl;

            for (int i=0; i<25; i++)
            {
                cout<< setw(12) << rnd()
                    << setw(12) << long(::lrand48())
                    << setw(12) << long(mtwist())
                    << setw(12) << long(deflt()) << nl;
            }
        }
    }

    {
        Random rnd(123456);
        Rand48 manual(123456);
        mtwist.seed(123456);
        deflt.seed(123456);

        // Two passes to ensure that reset is working properly
        for (label iter=0; iter < maxIter; ++iter)
        {
            ::srand48(123456);
            rnd.reset(123456);
            manual.seed(123456);
            mtwist.seed(123456);
            deflt.seed(123456);

            cout<< nl << "Random (Rand48) with seed = " << rnd.seed()
                << " interval [0,1000]" << nl;

            cout<< setw(12) << "Rand48()"
                << setw(12) << "drand48()";

            #ifdef TEST_RAW_IEEE
            cout<< setw(12) << "manual";
            #endif
            cout<< setw(12) << "mtwister";
            cout<< nl;

            for (int i=0; i<25; i++)
            {
                cout<< setw(12) << (rnd.sample01<scalar>()*1000)
                    << setw(12) << (drand48()*1000);

                #ifdef TEST_RAW_IEEE
                cout<< setw(12) << (randomFraction(manual.raw())*1000);
                #endif
                cout<< setw(12) << (uniform01(mtwist)*1000);
                cout<< setw(12) << (uniform01(deflt)*1000);
                cout<< nl;
            }
        }
    }

    {
        Rand48 rnd1(123456);
        Rand48 rnd2(123456);
        ::srand48(123456);

        rnd2.discard(10);

        cout<< nl << "Rand48 - test with offset of 10" << nl;

        cout<< setw(12) << "Rand48()"
            << setw(12) << "offset-10"
            << setw(12) << "lrand48()"
            << nl;

        for (int i=0; i<25; i++)
        {
            cout<< setw(12) << (rnd1())
                << setw(12) << (rnd2())
                << setw(12) << long(::lrand48())
                << nl;
        }
    }

    // Test uniformity of random
    testPosition(20);
    testPosition(3);

    // This should fail (in FULLDEBUG)
    const bool throwingError = FatalError.throwExceptions();
    try
    {
        Info<<"Random position(10,5): "
            << Random().position<label>(10, 5) << endl;
    }
    catch (Foam::error& err)
    {
        Info<< "Caught FatalError " << err << nl << endl;
    }

    FatalError.throwExceptions(throwingError);

    Info<< "\nDone" << nl << endl;

    return 0;
}


// ************************************************************************* //
