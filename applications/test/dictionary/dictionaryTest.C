/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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
    dictionaryTest

Description

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "IFstream.H"
#include "dictionary.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    IFstream dictStream("testDict");
    dictionary testDict(dictStream);

    Info<< testDict << endl;

    {
        dictionary someDict;
        someDict.add(keyType("a.*", true), "subdictValue");

        dictionary dict;
        dict.add("someDict", someDict);
        dict.add(keyType(".*", true), "parentValue");

        Info<< "dict:" << dict << endl;

        // Wildcard find.
        Info<< "Wildcard find \"abc\" in top directory : "
            << dict.lookup("abc") << endl;
        Info<< "Wildcard find \"abc\" in sub directory : "
            << dict.subDict("someDict").lookup("abc")
            << endl;
        Info<< "Recursive wildcard find \"def\" in sub directory : "
            << dict.subDict("someDict").lookup("def", true)
            << endl;
    }

    return 0;
}


// ************************************************************************* //
