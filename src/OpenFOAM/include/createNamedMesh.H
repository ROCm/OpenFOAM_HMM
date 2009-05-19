    word regionName;

    if (args.optionReadIfPresent("region", regionName))
    {
        Info<< "Create mesh " << regionName << " for time = "
            << runTime.timeName() << nl << endl;
    }
    else
    {
        regionName = fvMesh::defaultRegion;
        Info<< "Create mesh for time = "
            << runTime.timeName() << nl << endl;
    }

    fvMesh mesh
    (
        IOobject
        (
            regionName,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );
