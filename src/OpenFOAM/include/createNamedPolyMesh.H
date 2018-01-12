Foam::word regionName = Foam::polyMesh::defaultRegion;

if (args.readIfPresent("region", regionName))
{
    Foam::Info
        << "Create polyMesh " << regionName << " for time = "
        << runTime.timeName() << Foam::nl << Foam::endl;
}
else
{
    Foam::Info
        << "Create polyMesh for time = "
        << runTime.timeName() << Foam::nl << Foam::endl;
}

Foam::polyMesh mesh
(
    Foam::IOobject
    (
        regionName,
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    )
);
