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
    Test-dictionary2

Description

    Test dictionary insertion

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOstreams.H"
#include "IOobject.H"
#include "IFstream.H"
#include "dictionary.H"
#include "stringOps.H"

using namespace Foam;

void entryInfo(entry* e)
{
    if (e)
    {
        Info<<"added "
            << e->keyword() << ": " << typeid(e).name() << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();

    argList args(argc, argv);

    dictionary dict1;
    for (label i=0; i<5; ++i)
    {
        dictionary tmpdict;

        {
            entry* e = dict1.add
            (
                Foam::name("entry%d", i),
                string("entry" + Foam::name(i))
            );
            entryInfo(e);
        }

        {
            entry* e = tmpdict.add
            (
                Foam::name("subentry%d", i),
                string("subentry" + Foam::name(i))
            );
            entryInfo(e);
        }

        {
            entry* e = dict1.add
            (
                Foam::name("dict%d", i),
                tmpdict
            );
            entryInfo(e);
        }
    }

    // Insert new dictionary or merge into existing one
    for (auto k : { "dict1", "dict10" })
    {
        const word key(k);
        entry* e = dict1.add(key, dictionary(), true);

        if (e && e->isDict())
        {
            e->dict().add(word("sub1" + key), 10);
            e->dict().add(word("sub2" + key), 20);
            e->dict().add(word("sub3" + key), 30);
            e->dict().add(word("sub4" + key), 40);
            e->dict().add(word("sub5" + key), 50);
            e->dict().add(word("sub6" + key), 60);
        }
    }


    // overwrite existing
    {
        dict1.set("entry3", 1000); // overwrite
        entry* e = dict1.set(word("dict3"), 1000); // overwrite
        entryInfo(e);
    }

    // merge into existing dictionary: returns pointer to existing dict
    {
        dictionary tmpdict;
        tmpdict.add(word("something"), 3.14159);

        entry* e = dict1.add(word("dict4"), tmpdict, true);  // merge
        entryInfo(e);

        if (e) Info<< nl << "=> " << *e << nl;

        tmpdict.clear();
        tmpdict.add(word("other"), 2.718281);

        dict1.add(word("dict1"), tmpdict, true);  // merge
    }

    Info<< nl << "dictionary" << nl << nl;
    dict1.write(Info, false);


    {
        dict1.foundCompat
        (
            "newEntry", {{"entry1", 1612}, {"entry15", 1606}}
        );
        dict1.foundCompat
        (
            "newEntry", {{"entry15", 1612}, {"entry2", 1606}}
        );
        dict1.foundCompat
        (
            "newEntry", {{"entry3", 240}, {"entry2", 1606}}
        );

        // And some success
        dict1.foundCompat
        (
            "entry4", {{"none", 240}, {"entry2", 1606}}
        );
    }

    {
        label lval = readLabel
        (
            dict1.lookupCompat
            (
                "entry400", {{"none", 240}, {"entry3", 1606}}
            )
        );
        Info<< "int value: " << lval << nl;
    }


    // Could have different dictionary names and different entries names etc.
    // Quite ugly!
    {
        scalar sval = readScalar
        (
            dict1.csearchCompat
            (
                "newdictName", {{"dict4", 1706}, {"dict1", 1606}}
            )
            .dict()
            .lookupCompat
            (
                "newval", {{"something", 1606}, {"other", 1612}}
            )
        );

        Info<< "scalar value: " << sval << nl;

        sval = readScalar
        (
            dict1.csearchCompat
            (
                "newdictName", {{"dict1", 1606}, {"dict4", 1706}}
            )
            .dict()
            .lookupCompat
            (
                "newval", {{"something", 1606}, {"other", 1612}}
            )
        );

        Info<< "scalar value = " << sval << nl;
    }


    Info<< nl << "dictionary" << nl << nl;
    dict1.write(Info, false);

    Info<< "\nDone\n" << endl;

    return 0;
}


// ************************************************************************* //
