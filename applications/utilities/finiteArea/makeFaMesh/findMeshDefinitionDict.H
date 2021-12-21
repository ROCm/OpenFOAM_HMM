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
    Search for the appropriate faMeshDefinition dictionary...

\*---------------------------------------------------------------------------*/

const word dictName("faMeshDefinition");

autoPtr<IOdictionary> meshDictPtr;

{
    fileName dictPath;

    const word& regionDir =
        (regionName == polyMesh::defaultRegion ? word::null : regionName);

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
        // Check global location
        exists
        (
            runTime.path()/runTime.caseConstant()
            /regionDir/faMesh::meshSubDir/dictName
        )
    )
    {
        // Dictionary present in constant faMesh directory (old-style)

        dictPath =
            runTime.constant()
            /regionDir/faMesh::meshSubDir/dictName;

        // Warn that constant/faMesh/faMeshDefinition was used
        // instead of system/faMeshDefinition
        #if 0
        WarningIn(args.executable())
            << "Using the old faMeshDefinition location: "
            << dictPath << nl
            << "    instead of default location: "
            << runTime.system()/regionDir/dictName << nl
            << endl;
        #endif
    }
    else
    {
        // Assume dictionary is in the system directory

        dictPath = runTime.system()/regionDir/dictName;
    }

    IOobject meshDictIO
    (
        dictPath,
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false,  // no registerObject
        true    // is globalObject
    );

    if (!meshDictIO.typeHeaderOk<IOdictionary>(true))
    {
        FatalErrorInFunction
            << meshDictIO.objectPath() << nl
            << exit(FatalError);
    }

    Info<< "Creating faMesh from definition: "
        << meshDictIO.objectRelPath() << endl;

    meshDictPtr = autoPtr<IOdictionary>::New(meshDictIO);
}

IOdictionary& meshDefDict = *meshDictPtr;


// ************************************************************************* //
