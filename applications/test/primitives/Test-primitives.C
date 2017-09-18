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

Application
    Test-primitives

Description
    Parsing etc for primitives.

\*---------------------------------------------------------------------------*/

#include "scalar.H"
#include "label.H"
#include "StringStream.H"
#include "NASCore.H"
#include "parsing.H"
#include "Tuple2.H"

using namespace Foam;

// Shadow fileFormats::NASCore::readNasScalar
inline scalar readNasScalar(const std::string& str)
{
    return fileFormats::NASCore::readNasScalar(str);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class TYPE>
unsigned testParsing
(
    TYPE (*function)(const std::string&),
    const List<Tuple2<std::string, bool>>& tests
)
{
    unsigned nFail = 0;

    // Expect some failures
    const bool prev = FatalIOError.throwExceptions();

    for (const Tuple2<std::string, bool>& test : tests)
    {
        const std::string& str = test.first();
        const bool expected = test.second();

        bool parsed = true;

        TYPE val;
        try
        {
            val = function (str);
        }
        catch (Foam::error& err)
        {
            parsed = false;
        }

        if (parsed)
        {
            if (expected)
            {
                Info<< "(pass) parsed " << str << " = " << val << nl;
            }
            else
            {
                ++nFail;
                Info<< "(fail) unexpected success for " << str << nl;
            }
        }
        else
        {
            if (expected)
            {
                ++nFail;
                Info<< "(fail) unexpected failure " << str << nl;
            }
            else
            {
                Info<< "(pass) expected failure " << str << nl;
            }
        }
    }

    FatalIOError.throwExceptions(prev);

    return nFail;
}


int main(int argc, char *argv[])
{
    unsigned nFail = 0;

    {
        Info<< nl << "Test readDouble:" << nl;
        nFail += testParsing
        (
            &readDouble,
            {
                { "", false },
                { "  ", false },
                { " xxx ", false },
                { " 1234E-", false },
                { " 1234E junk", false },
                { " 3.14159 ", true },
                { " 31.4159E-1 " , true },
            }
        );
    }

    {
        Info<< nl << "Test readFloat:" << nl;
        nFail += testParsing
        (
            &readFloat,
            {
                { " 3.14159 ", true },
                { " 31.4159E-1 " , true },
                { " 31.4159E200 " , false },
                { " 31.4159E20 " , true },
            }
        );
    }

    {
        Info<< nl << "Test readNasScalar:" << nl;
        nFail += testParsing
        (
            &readNasScalar,
            {
                { " 3.14159 ", true },
                { " 31.4159E-1 " , true },
                { " 314.159-2 " , true },
                { " 31.4159E200 " , true },
                { " 31.4159E20 " , true },
            }
        );
    }

    {
        Info<< nl << "Test readInt32 (max= " << INT32_MAX << "):" << nl;
        nFail += testParsing
        (
            &readInt32,
            {
                { " 3.14159 ", false },
                { " 31.4159E-1 " , false },
                { "100" , true },
                { "	2147483644" , true },
                { "   2147483700  " , false },
            }
        );
    }

    {
        Info<< nl << "Test readUint32 (max= " << INT32_MAX << "):" << nl;
        nFail += testParsing
        (
            &readUint32,
            {
                { "	2147483644" , true },
                { "   2147483700  " , true },
            }
        );
    }

    if (nFail)
    {
        Info<< nl << "failed " << nFail << " tests" << nl;
        return 1;
    }

    Info<< nl << "passed all tests" << nl;
    return 0;
}

// ************************************************************************* //
