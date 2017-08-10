/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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
    Test-faces

Description
    Simple tests for various faces

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "labelledTri.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    face f1{ 1, 2, 3, 4 };
    Info<< "face:" << f1 << nl;

    triFace t1{ 1, 2, 3 };
    Info<< "triFace:" << t1 << nl;

    f1 = t1;
    Info<< "face:" << f1 << nl;

    f1 = t1.triFaceFace();
    Info<< "face:" << f1 << nl;

    // expect these to fail
    const bool throwingError = FatalError.throwExceptions();
    try
    {
        labelledTri l1{ 1, 2, 3, 10, 24 };
        Info<< "labelled:" << l1 << nl;
    }
    catch (Foam::error& err)
    {
        WarningInFunction
            << "Caught FatalError " << err << nl << endl;
    }
    FatalError.throwExceptions(throwingError);

    labelledTri l2{ 1, 2, 3 };
    Info<< "labelled:" << l2 << nl;

    labelledTri l3{ 1, 2, 3, 10 };
    Info<< "labelled:" << l3 << nl;

    t1.flip();
    l3.flip();

    Info<< "flip:" << t1 << nl;
    Info<< "flip:" << l3 << nl;

    return 0;
}


// ************************************************************************* //
