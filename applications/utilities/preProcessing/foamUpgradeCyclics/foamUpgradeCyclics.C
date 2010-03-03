/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    foamUpgradeCyclics

Description
    Simple tool to upgrade mesh and fields for split cyclics

Usage

    - foamUpgradeCyclics [OPTION]

    @param -test \n
    Suppress writing the updated files with split cyclics

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IOdictionary.H"
#include "polyMesh.H"
#include "entry.H"
#include "IOPtrList.H"
#include "cyclicPolyPatch.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<entry>, 0);
}

// Main program:

int main(int argc, char *argv[])
{
    argList::addBoolOption("test");
    #include "addRegionOption.H"

#   include "setRootCase.H"
#   include "createTime.H"

    Foam::word regionName = polyMesh::defaultRegion;
    args.optionReadIfPresent("region", regionName);

    fileName regionPrefix = "";
    if (regionName != polyMesh::defaultRegion)
    {
        regionPrefix = regionName;
    }


    // Per cyclic patch the new name for this side and the other side
    HashTable<word> thisNames;
    HashTable<word> nbrNames;

    // Read boundary file without reading mesh
    {
        IOobject io
        (
            "boundary",
            runTime.findInstance
            (
                regionPrefix/polyMesh::meshSubDir,
                "boundary"
            ),
            polyMesh::meshSubDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );
        

        Info<< "Reading boundary from " << io.filePath() << endl;


        // Read PtrList of dictionary.
        const word oldTypeName = IOPtrList<entry>::typeName;
        const_cast<word&>(IOPtrList<entry>::typeName) = word::null;
        IOPtrList<entry> patches(io);
        const_cast<word&>(IOPtrList<entry>::typeName) = oldTypeName;
        // Fake type back to what was in field
        const_cast<word&>(patches.type()) = patches.headerClassName();

        // Temporary convert to dictionary
        dictionary patchDict;
        forAll(patches, i)
        {
            patchDict.add(patches[i].keyword(), patches[i].dict());
        }

        // Replace any 'cyclic'
        label nOldCyclics = 0;
        forAll(patches, patchI)
        {
            const dictionary& patchDict = patches[patchI].dict();

            if (word(patchDict["type"]) == cyclicPolyPatch::typeName)
            {
                if (patchDict.found("neighbourPatch"))
                {
                    Info<< "Patch " << patches[patchI].keyword()
                        << " already has 'neighbourPatch' entry; assuming it"
                        << " is already converted." << endl;
                }
                else
                {
                    Info<< "Patch " << patches[patchI].keyword()
                        << " does not have 'neighbourPatch' entry; assuming it"
                        << " is of the old type." << endl;
                    nOldCyclics++;
                }
            }
        }

        Info<< "Detected " << nOldCyclics << " old cyclics." << nl << endl;


        // edo the 



        PtrList<entry> oldPatches(patches);

Pout<< "oldPatches:" << oldPatches << endl;

        // Extend
        label nOldPatches = patches.size();
        patches.setSize(nOldPatches+nOldCyclics);


        // Add new entries
        label newPatchI = 0;
        forAll(oldPatches, patchI)
        {
            const dictionary& patchDict = oldPatches[patchI].dict();

            if
            (
                word(patchDict["type"]) == cyclicPolyPatch::typeName
            && !patchDict.found("neighbourPatch")
            )
            {
                const word& name = oldPatches[patchI].keyword();
                label nFaces = readLabel(patchDict["nFaces"]);
                label startFace = readLabel(patchDict["startFace"]);

                word thisName = name + "_half0";
                word nbrName = name + "_half1";

                thisNames.insert(name, thisName);
                nbrNames.insert(name, nbrName);

                // Change entry on this side
                patches.set(newPatchI, oldPatches(patchI));
                dictionary& thisPatchDict = patches[newPatchI].dict();
                thisPatchDict.add("neighbourPatch", nbrName);
                thisPatchDict.set("nFaces", nFaces/2);
                patches[newPatchI].keyword() = thisName;
                newPatchI++;

                // Add entry on other side
                patches.set(newPatchI, oldPatches(patchI));
                dictionary& nbrPatchDict = patches[newPatchI].dict();
                nbrPatchDict.add("neighbourPatch", nbrName);
                nbrPatchDict.set("nFaces", nFaces/2);
                nbrPatchDict.set("startFace", startFace+nFaces/2);
                patches[newPatchI].keyword() = nbrName;
                newPatchI++;
            }
            else
            {
                patches.set(newPatchI++, oldPatches(patchI));
            }
        }

        Info<< "boundary:" << patches << endl;

        if (returnReduce(nOldCyclics, sumOp<label>()) >= 0)
        {
            if (args.optionFound("test"))
            {
                Info<< "-test option: no changes made" << nl << endl;
            }
            else
            {
                if (mvBak(patches.objectPath(), "old"))
                {
                    Info<< "Backup to    "
                        << (patches.objectPath() + ".old") << nl;
                }

                Info<< "Write  to    "
                    << patches.objectPath() << nl << endl;
                patches.write();
            }
        }
    }


//    {
//        // Read dictionary. (disable class type checking so we can load
//        // field)
//        Info<< "Loading dictionary " << fieldName << endl;
//        const word oldTypeName = IOdictionary::typeName;
//        const_cast<word&>(IOdictionary::typeName) = word::null;
//
//        IOdictionary fieldDict
//        (
//            IOobject
//            (
//                "p",
//                instance,
//                mesh,
//                IOobject::MUST_READ,
//                IOobject::NO_WRITE,
//                false
//            )
//        );
//        const_cast<word&>(IOdictionary::typeName) = oldTypeName;
//        // Fake type back to what was in field
//        const_cast<word&>(fieldDict.type()) = fieldDict.headerClassName();
//
//        Info<< "Loaded dictionary " << fieldName
//            << " with entries " << fieldDict.toc() << endl;
//
//        dictionary& boundaryField = fieldDict.subDict("boundaryField");
//
//        forAllConstIter(HashTable<word>, thisNames, iter)
//        {
//            const word& patchName = iter.key();
//            const word& newName = iter();
//
//            Info<< "Looking for entry for patch " << patchName << endl;
//
//            if (boundaryField.found(patchName) && !boundaryField.found(iter()))
//            {
//                const dictionary& patchDict = boundaryField[patchName];
//
//                Field<scalar> fld(patchDict.lookup("value"));
//
//                
//            }
//
//
//        forAllIter(IDLList<entry>, boundaryField, patchIter)
//        {
//            
//        }


    return 0;
}


// ************************************************************************* //
