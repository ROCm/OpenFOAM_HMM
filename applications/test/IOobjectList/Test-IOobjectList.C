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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addOption("re", "wordReList");

    // timeSelector::addOptions();
    timeSelector::addOptions(true, true);

    #include "setRootCase.H"
    #include "createTime.H"

    wordReList matcher;
    if (args.optionFound("re"))
    {
        matcher = args.optionReadList<wordRe>("re");
        Info<<"limit names: " << matcher << nl;

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

        Info<<"Name:    " << flatOutput(objects.sortedNames()) << nl
            <<"Objects: " << objects << nl
            <<"Classes: " << classes << nl;

        classes.filterKeys(subsetTypes);
        Info<<"only retain: " << flatOutput(subsetTypes) << nl;
        Info<<"Pruned: " << classes << nl;

        classes = objects.classes();
        classes.erase(subsetTypes);
        Info<<"remove: " << flatOutput(subsetTypes) << nl;
        Info<<"Pruned: " << classes << nl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
