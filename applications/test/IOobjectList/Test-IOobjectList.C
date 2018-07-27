/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenCFD Ltd.
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
    Basic tests of IOobjectList

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "volFields.H"
#include "timeSelector.H"
#include "IOobjectList.H"
#include "hashedWordList.H"

using namespace Foam;

void report(const IOobjectList& objects)
{
    Info<< "Names: " << flatOutput(objects.sortedNames()) << nl
        << "Objects: " << objects << nl
        << "----" << nl;
}


void reportDetail(const IOobjectList& objects)
{
    Info<<"Details:" << nl;

    for (const word& key : objects.sortedNames())
    {
        IOobject* io = objects.lookup(key);

        Info<< key << " (" << io->headerClassName()
            << ") = addr " << long(io) << nl;
    }

    Info<<"====" << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addOption
    (
        "filter",
        "wordRes",
        "filter keys with names or regexs"
    );
    argList::addBoolOption
    (
        "copy-append",
        "test move append lists (requires -filter)"
    );
    argList::addBoolOption
    (
        "move-append",
        "test move append lists (requires -filter)"
    );

    // timeSelector::addOptions();
    timeSelector::addOptions(true, true);

    #include "setRootCase.H"
    #include "createTime.H"

    wordRes matcher;
    if (args.readListIfPresent<wordRe>("filter", matcher))
    {
        Info<<"limit names: " << matcher << nl;
    }

    if (args.found("copy-append") && matcher.empty())
    {
        FatalError
            << nl << "The -copy-append test also requires -filter" << nl
            << exit(FatalError);
    }
    if (args.found("move-append") && matcher.empty())
    {
        FatalError
            << nl << "The -move-append test also requires -filter" << nl
            << exit(FatalError);
    }


    const hashedWordList subsetTypes
    {
        volScalarField::typeName,
        volScalarField::Internal::typeName,
        volVectorField::typeName,
    };


    instantList timeDirs = timeSelector::select0(runTime, args);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        // Objects at this time
        IOobjectList objects(runTime, runTime.timeName());
        HashTable<wordHashSet> classes =
        (
            matcher.size()
          ? objects.classes(matcher)
          : objects.classes()
        );

        Info<< "Time: " << runTime.timeName() << nl;

        report(objects);

        classes.filterKeys(subsetTypes);
        Info<<"only retain: " << flatOutput(subsetTypes) << nl;
        Info<<"Pruned: " << classes << nl;

        classes = objects.classes();
        classes.erase(subsetTypes);
        Info<<"remove: " << flatOutput(subsetTypes) << nl;
        Info<<"Pruned: " << classes << nl;

        // On last time
        if (timeI == timeDirs.size()-1)
        {
            if (args.found("copy-append"))
            {
                Info<< nl << "Test move append" << nl;
            }
            else if (args.found("move-append"))
            {
                Info<< nl << "Test move append" << nl;
            }
            else
            {
                continue;
            }

            IOobjectList other(runTime, runTime.timeName());

            Info<< "==original==" << nl; reportDetail(objects);

            objects.filterKeys(matcher);

            Info<< "==target==" << nl; reportDetail(objects);
            Info<< "==source==" << nl; reportDetail(other);

            if (args.found("copy-append"))
            {
                objects.append(other);

                Info<< nl << "After copy-append" << nl;
            }
            else
            {
                objects.append(std::move(other));

                Info<< nl << "After move-append" << nl;
            }

            Info<< "==target==" << nl; reportDetail(objects);
            Info<< "==source==" << nl; reportDetail(other);

            Info<< nl;
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
