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
    std::initializer_list
    <
        Tuple2<bool, std::string>
    > tests
)
{
    unsigned nFail = 0;
    string errMsg;

    // Expect some failures
    const bool prev = FatalIOError.throwExceptions();

    for (const Tuple2<bool, std::string>& test : tests)
    {
        const bool expected = test.first();
        const std::string& str = test.second();

        bool parsed = true;

        TYPE val;
        try
        {
            val = function (str);
        }
        catch (Foam::error& err)
        {
            parsed = false;
            errMsg = err.message();
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
                Info<< "(fail) unexpected";
            }
            else
            {
                Info<< "(pass) expected";
            }

            Info<< " failure " << str
                << "  >> " << errMsg.c_str() << nl;
        }
    }

    FatalIOError.throwExceptions(prev);

    return nFail;
}


int main(int argc, char *argv[])
{
    unsigned nFail = 0;

    {
        Info<< nl << "Test readDouble: (small=" << doubleScalarVSMALL
            << " great=" << doubleScalarVSMALL << "):" << nl;
        nFail += testParsing
        (
            &readDouble,
            {
                { false, "" },
                { false, "  " },
                { false, " xxx " },
                { false, " 1234E-" },
                { false, " 1234E junk" },
                { true,  " 3.14159 " },
                { true,  " 31.4159E-1 "  },
                { false, " 100E1000 "  },
                { true,  " 1E-40 "  },
                { true,  " 1E-305 "  },
                { true,  " 1E-37 "  },
                { true,  " 1E-300 "  },
            }
        );
    }

    {
        Info<< nl << "Test readFloat: (small=" << floatScalarVSMALL
            << " great=" << floatScalarVGREAT << "):" << nl;

        nFail += testParsing
        (
            &readFloat,
            {
                { true,  " 3.14159 " },
                { true,  " 31.4159E-1 "  },
                { false, " 31.4159E200 "  },
                { true,  " 31.4159E20 "  },
                { true,  " 1E-40 "  },
                { true,  " 1E-305 "  },
                { true,  " 1E-37 "  },
                { true,  " 1E-300 "  },
            }
        );
    }

    {
        Info<< nl << "Test readNasScalar:" << nl;
        nFail += testParsing
        (
            &readNasScalar,
            {
                { true,  " 3.14159 " },
                { true,  " 31.4159E-1 "  },
                { true,  " 314.159-2 "  },
                { true,  " 31.4159E200 "  },
                { true,  " 31.4159E20 "  },
                { true,  " 1E-40 "  },
                { true,  " 1E-305 "  },
                { true,  " 1E-37 "  },
                { true,  " 1E-300 "  },
            }
        );
    }

    {
        Info<< nl << "Test readInt32 (max=" << INT32_MAX << "):" << nl;
        nFail += testParsing
        (
            &readInt32,
            {
                { false, " 3.14159 " },
                { false, " 31E1 " },
                { false, " 31.4159E-1 "  },
                { true,  "100"  },
                { true,  "	2147483644"  },
                { false, "   2147483700  "  },
            }
        );
    }

    {
        Info<< nl << "Test readUint32 (max="
            << unsigned(UINT32_MAX) << "):" << nl;
        nFail += testParsing
        (
            &readUint32,
            {
                { true,  "\t2147483644"  },
                { true,  " 2147483700  "  },
                { true,  " 4294967295  "  },
                { false, " 4294968000  "  },
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
