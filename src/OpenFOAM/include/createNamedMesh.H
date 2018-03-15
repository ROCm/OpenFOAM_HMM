Foam::word regionName = Foam::fvMesh::defaultRegion;

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

Foam::fvMesh mesh
(
    Foam::IOobject
    (
        regionName,
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    )
);
