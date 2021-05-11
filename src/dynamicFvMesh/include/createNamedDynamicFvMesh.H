    Foam::word regionName = Foam::polyMesh::defaultRegion;

    if (args.readIfPresent("region", regionName))
    {
        Foam::Info
            << "Create mesh " << regionName << " for time = "
            << runTime.timeName() << Foam::nl << Foam::endl;
    }
    else
    {
        Foam::Info
            << "Create mesh for time = "
            << runTime.timeName() << Foam::nl << Foam::endl;
    }

    autoPtr<dynamicFvMesh> meshPtr
    (
        dynamicFvMesh::New
        (
            IOobject
            (
                regionName,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );

    dynamicFvMesh& mesh = meshPtr();
