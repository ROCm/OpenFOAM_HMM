/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Create a dynamicFvMesh for a specified region, or the defaultRegion

Required Variables
    - runTime [Time]

Provided Variables
    - regionName [word]
    - mesh [dynamicFvMesh], meshPtr

\*---------------------------------------------------------------------------*/

Foam::word regionName(Foam::polyMesh::defaultRegion);
Foam::autoPtr<Foam::dynamicFvMesh> meshPtr;

{
    Foam::Info << "Create dynamic mesh";
    if (args.readIfPresent("region", regionName))
    {
        Foam::Info << ' ' << regionName;
    }
    Foam::Info << " for time = " << runTime.timeName() << Foam::nl;

    meshPtr = dynamicFvMesh::New
    (
        IOobject
        (
            regionName,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );
}


dynamicFvMesh& mesh = meshPtr();

Foam::Info << Foam::endl;


// ************************************************************************* //
