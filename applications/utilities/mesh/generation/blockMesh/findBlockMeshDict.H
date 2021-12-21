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
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Search for the appropriate blockMeshDict dictionary...

\*---------------------------------------------------------------------------*/

// Search for the appropriate blockMesh dictionary....
const word dictName("blockMeshDict");

autoPtr<IOdictionary> meshDictPtr;

{
    fileName dictPath;

    if (args.readIfPresent("dict", dictPath))
    {
        // Dictionary specified on the command-line ...

        if (isDir(dictPath))
        {
            dictPath /= dictName;
        }
    }
    else if
    (
        exists
        (
            runTime.path()/runTime.constant()
            /regionPath/polyMesh::meshSubDir/dictName
        )
    )
    {
        // Dictionary present in constant polyMesh directory (old-style)

        dictPath =
            runTime.constant()
                /regionPath/polyMesh::meshSubDir/dictName;


        // Warn that constant/polyMesh/blockMeshDict was used
        // instead of system/blockMeshDict
        WarningIn(args.executable())
            << "Using the old blockMeshDict location: "
            << dictPath << nl
            << "    instead of the default location:  "
            << runTime.system()/regionPath/dictName << nl
            << endl;
    }
    else
    {
        // Assume dictionary is to be found in the system directory

        dictPath = runTime.system()/regionPath/dictName;
    }

    IOobject meshDictIO
    (
        dictPath,
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (!meshDictIO.typeHeaderOk<IOdictionary>(true))
    {
        FatalErrorInFunction
            << meshDictIO.objectPath() << nl
            << exit(FatalError);
    }

    Info<< "Creating block mesh from "
        << meshDictIO.objectRelPath() << endl;

    meshDictPtr = autoPtr<IOdictionary>::New(meshDictIO);
}

const IOdictionary& meshDict = *meshDictPtr;


// ************************************************************************* //
