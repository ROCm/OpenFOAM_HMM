/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Description

\*---------------------------------------------------------------------------*/

#include "NamedEnum.H"
#include "IOstreams.H"

using namespace Foam;

class namedEnumTest
{
public:

    enum option
    {
        a,
        b,
        c,
        d
    };

    static const Foam::NamedEnum<option, 4> namedEnum;
};


template<>
const char* Foam::NamedEnum<namedEnumTest::option, 4>::names[] =
{
    "a",
    "b",
    "c",
    "d"
};

const Foam::NamedEnum<namedEnumTest::option, 4> namedEnumTest::namedEnum;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    const List<namedEnumTest::option> options
        = namedEnumTest::namedEnum.enums();

    Info<< "enums: " << options << nl;

    Info<< "loop over enums (as list):" << nl;
    forAll(options, i)
    {
        const namedEnumTest::option& opt = options[i];

        Info<< "option[" << opt
            << "] = '" << namedEnumTest::namedEnum[opt] << "'" << nl;
    }

#if __cplusplus > 201100L
    // C++11
    Info<< "loop over enums (C++11 for range):" << nl;
    for (auto const& opt : options)
    {
        Info<< "option[" << opt
            << "] = '" << namedEnumTest::namedEnum[opt] << "'" << nl;
    }
#else
    Info<< "loop over enums (via iterator):" << nl;
    forAllConstIter(List<namedEnumTest::option>, options, iter)
    {
        const namedEnumTest::option& opt = *iter;

        Info<< "option[" << opt
            << "] = '" << namedEnumTest::namedEnum[opt] << "'" << nl;
    }
#endif

    Info<< nl
        << namedEnumTest::namedEnum["a"] << nl
        << namedEnumTest::namedEnum[namedEnumTest::a] << nl;

    Info<< "--- test read construction ---" << endl;

    namedEnumTest::option dummy(namedEnumTest::namedEnum.read(Sin));
    Info<< namedEnumTest::namedEnum[dummy] << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
