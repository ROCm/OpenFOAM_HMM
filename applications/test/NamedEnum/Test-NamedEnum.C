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

    enum class option
    {
        A,
        B,
        C,
        D
    };

    static const Foam::NamedEnum<option, 4> optionNamed;
};


template<>
const char* Foam::NamedEnum<namedEnumTest::option, 4>::names[] =
{
    "a",
    "b",
    "c",
    "d",
};

const Foam::NamedEnum<namedEnumTest::option, 4> namedEnumTest::optionNamed;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    const List<namedEnumTest::option> options
        = namedEnumTest::namedEnum.enums();

    dictionary testDict;
    testDict.add("lookup1", "c");

    Info<< nl
        << int(namedEnumTest::optionNamed["a"]) << nl
        << namedEnumTest::optionNamed[namedEnumTest::option::A] << nl;

    Info<< "--- test dictionary lookup ---" << endl;
    {
        Info<< "dict: " << testDict << endl;

        Info<< "got: "
            <<  int
                (
                    namedEnumTest::optionNamed.lookupOrDefault
                    (
                        "notFound",
                        testDict,
                        namedEnumTest::option::A
                    )
                )
            << nl;

        Info<< "got: "
            <<  int
                (
                    namedEnumTest::optionNamed.lookupOrDefault
                    (
                        "lookup1",
                        testDict,
                        namedEnumTest::option::A
                    )
                )
            << nl;
    }

    Info<< "--- test read ---" << endl;

    namedEnumTest::option dummy(namedEnumTest::optionNamed.read(Sin));
    Info<< namedEnumTest::optionNamed[dummy] << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
