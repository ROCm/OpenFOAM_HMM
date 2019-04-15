Foam::autoPtr<Foam::fvMesh> meshPtr(nullptr);
Foam::word regionName = Foam::fvMesh::defaultRegion;

if (args.found("dry-run") || args.found("dry-run-write"))
{
    Foam::Info
        << "Operating in 'dry-run' mode: case will run for 1 time step.  "
        << "All checks assumed OK on a clean exit" << Foam::endl;

    // Allow region in combination with dry-run
    args.readIfPresent("region", regionName);

    Foam::FieldBase::allowConstructFromLargerSize = true;

    // Create a simplified 1D mesh and attempt to re-create boundary conditions
    meshPtr.reset
    (
        new Foam::simplifiedMeshes::columnFvMesh(runTime, regionName)
    );

    // Stop after 1 iteration of the simplified mesh

    if (args.found("dry-run-write"))
    {
        // Using saWriteNow triggers function objects execute(), write()
        runTime.stopAt(Foam::Time::saWriteNow);

        // Make sure mesh gets output to the current time (since instance
        // no longer constant)
        meshPtr().setInstance(runTime.timeName());
    }
    else
    {
        // Using saNoWriteNow triggers function objects execute(),
        // but not write()
        runTime.stopAt(Foam::Time::saNoWriteNow);
    }

    Foam::functionObject::outputPrefix = "postProcessing-dry-run";
}
else
{
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

    meshPtr.reset
    (
        new Foam::fvMesh
        (
            Foam::IOobject
            (
                regionName,
                runTime.timeName(),
                runTime,
                Foam::IOobject::MUST_READ
            )
        )
    );
}

Foam::fvMesh& mesh = meshPtr();
