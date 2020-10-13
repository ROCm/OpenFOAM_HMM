/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Test token construct assign etc.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOobject.H"
#include "IOstreams.H"
#include "IFstream.H"
#include "StringStream.H"
#include "cpuTime.H"
#include "labelList.H"
#include "DynamicList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();

    argList args(argc, argv, false, true);

    token tok1;
    Info<< "construct null: " << tok1.info() << endl;

    tok1 = double(3.14159);
    Info<< "assign double: " << tok1.info() << endl;

    token tok2(tok1);
    Info<< "copy construct: " << tok2.info() << endl;

    tok1 = word("this-word");
    Info<< "assign word: " << tok1.info() << endl;

    token tok3(tok1);
    Info<< "copy construct: " << tok3.info() << endl;
    Info<< "orig: " << tok1.info() << endl;

    token tok4(std::move(tok1));
    Info<< "move construct: " << tok4.info() << endl;
    Info<< "orig: " << tok1.info() << endl;

    tok3 = tok4;
    Info<< "assign token: " << tok3.info() << endl;
    Info<< "orig: " << tok4.info() << endl;

    //
    // Compound
    //

    {
        // This version is good

        token ctok1(new token::Compound<labelList>(identity(10)));

        Info<< "compound token: " << ctok1.info() << nl << ctok1 << endl;
    }

    {
        // This also works, but not actually using the autoPtr directly

        autoPtr<token::Compound<labelList>> ptr
        (
            new token::Compound<labelList>(identity(10, -9))
        );

        token ctok1(ptr.release());  // release() not get()!

        Info<< "compound token: " << ctok1.info() << nl << ctok1 << endl;
    }

    #if 0
    {
        // This version will segfault.
        // The implicit pointer cast from autoPtr to pointer wracks havoc

        autoPtr<token::Compound<labelList>> ptr
        (
            new token::Compound<labelList>(identity(10))
        );

        token ctok1(ptr);

        Info<< "compound token: " << ctok1.info() << nl << ctok1 << endl;
    }
    #endif

    return 0;
}

// ************************************************************************* //
