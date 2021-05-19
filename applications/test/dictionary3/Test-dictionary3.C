/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Test-dictionary3

Description
    Test expressions and re-expansions

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOstreams.H"
#include "dictionary.H"
#include "vector.H"
#include "StringStream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();

    {
        IStringStream is
        (
            "value   10;"
            "scalar1 $value;"
            "scalar2 -$value;"

            // Use #eval expansion entirely
            "vector1 ${{vector($value, -$value, $value)}};"
            "vector2 ($value -$value $value);"
        );

        dictionary dict(is);

        Info<< "input dictionary:" << dict << nl;

        Info<< "value: " << dict.get<scalar>("value") << nl;

        Info<< "scalar1: " << dict.get<scalar>("scalar1") << nl;
        Info<< "scalar2: " << dict.get<scalar>("scalar2") << nl;

        Info<< "vector1: " << dict.get<vector>("vector1") << nl;
        Info<< "vector2: " << dict.get<vector>("vector2") << nl;
    }

    return 0;
}


// ************************************************************************* //
