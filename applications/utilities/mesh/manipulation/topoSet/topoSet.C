/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
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
    Operates on cellSets/faceSets/pointSets through a dictionary.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "topoSetSource.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "dict",
        "file",
        "specify an alternative dictionary for the topoSet dictionary"
    );
#   include "addRegionOption.H"
    argList::addBoolOption
    (
        "noSync",
        "do not synchronise selection across coupled patches"
    );

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createNamedPolyMesh.H"

    const bool noSync = args.optionFound("noSync");

    const word dictName("topoSetDict");

    fileName dictPath = dictName;
    if (args.optionFound("dict"))
    {
        dictPath = args["dict"];
        if (isDir(dictPath))
        {
            dictPath = dictPath / dictName;
        }
    }

    Info<< "Reading " << dictName << "\n" << endl;

    IOdictionary topoSetDict
    (
        (
            args.optionFound("dict")
          ? IOobject
            (
                dictPath,
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
          : IOobject
            (
                dictName,
                runTime.system(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        ) 
    );


    // Read set construct info from dictionary
    PtrList<dictionary> patchSources(topoSetDict.lookup("actions"));

    forAll(patchSources, i)
    {
        const dictionary& dict = patchSources[i];

        const word setName(dict.lookup("name"));
        const word actionName(dict.lookup("action"));
        const word setType(dict.lookup("type"));


        topoSetSource::setAction action = topoSetSource::toAction(actionName);

        autoPtr<topoSet> currentSet;
        if
        (
            (action == topoSetSource::NEW)
         || (action == topoSetSource::CLEAR)
        )
        {
            currentSet = topoSet::New(setType, mesh, setName, 10000);
            Info<< "Created set " << setName << endl;
        }
        else if (action == topoSetSource::REMOVE)
        {
            //?
        }
        else
        {
            currentSet = topoSet::New
            (
                setType,
                mesh,
                setName,
                IOobject::MUST_READ
            );
            Info<< "Read set " << setName << " with size "
                << currentSet().size() << endl;
        }



        // Handle special actions (clear, invert) locally, rest through sources.
        switch (action)
        {
            case topoSetSource::NEW:
            case topoSetSource::ADD:
            case topoSetSource::DELETE:
            {
                Info<< "    Applying source " << word(dict.lookup("source"))
                    << endl;
                autoPtr<topoSetSource> source = topoSetSource::New
                (
                    dict.lookup("source"),
                    mesh,
                    dict.subDict("sourceInfo")
                );

                source().applyToSet(action, currentSet());
                // Synchronize for coupled patches.
                if (!noSync) currentSet().sync(mesh);
                currentSet().write();
            }
            break;

            case topoSetSource::SUBSET:
            {
                Info<< "    Applying source " << word(dict.lookup("source"))
                    << endl;
                autoPtr<topoSetSource> source = topoSetSource::New
                (
                    dict.lookup("source"),
                    mesh,
                    dict.subDict("sourceInfo")
                );

                // Backup current set.
                autoPtr<topoSet> oldSet
                (
                    topoSet::New
                    (
                        setType,
                        mesh,
                        currentSet().name() + "_old2",
                        currentSet()
                    )
                );

                currentSet().clear();
                source().applyToSet(topoSetSource::NEW, currentSet());

                // Combine new value of currentSet with old one.
                currentSet().subset(oldSet());
                // Synchronize for coupled patches.
                if (!noSync) currentSet().sync(mesh);
                currentSet().write();
            }
            break;

            case topoSetSource::CLEAR:
                Info<< "    Clearing set" << endl;
                currentSet().clear();
                currentSet().write();
            break;

            case topoSetSource::INVERT:
                Info<< "    Inverting set" << endl;
                currentSet().invert(currentSet().maxSize(mesh));
                currentSet().write();
            break;

            default:
                WarningIn(args.executable())
                    << "Unhandled action " << action << endl;
            break;
        }

        if (currentSet.valid())
        {
            Info<< "    Set " << currentSet().name()
                << " now size " << currentSet().size()
                << endl;
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
