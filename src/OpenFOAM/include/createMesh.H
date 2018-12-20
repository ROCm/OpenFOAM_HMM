Foam::autoPtr<Foam::fvMesh> meshPtr(nullptr);

if (args.found("dry-run") || args.found("dry-run-write"))
{
    Foam::Info
        << "Operating in 'dry-run' mode: case will run for 1 time step.  "
        << "All checks assumed OK on a clean exit" << Foam::endl;

    Foam::FieldBase::allowConstructFromLargerSize = true;

    // Create a simplified 1D mesh and attempt to re-create boundary conditions
    meshPtr.reset(new Foam::simplifiedMeshes::columnFvMesh(runTime));

    // Stop after 1 iteration of the simplified mesh

    if (args.found("dry-run-write"))
    {
        // Using saWriteNow triggers function objects execute(), write()
        runTime.stopAt(Foam::Time::saWriteNow);
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
    Foam::Info
        << "Create mesh for time = "
        << runTime.timeName() << Foam::nl << Foam::endl;

    meshPtr.reset
    (
        new Foam::fvMesh
        (
            Foam::IOobject
            (
                Foam::fvMesh::defaultRegion,
                runTime.timeName(),
                runTime,
                Foam::IOobject::MUST_READ
            )
        )
    );
}

Foam::fvMesh& mesh = meshPtr();
