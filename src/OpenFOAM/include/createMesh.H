Foam::autoPtr<Foam::fvMesh> meshPtr(nullptr);

if (args.optionFound("dry-run") || args.optionFound("dry-run-write"))
{
    Foam::Info
        << "Operating in 'dry-run' mode: case will run for 1 time step.  "
        << "All checks assumed OK on a clean exit" << Foam::endl;

    Foam::FieldBase::allowConstructFromLargerSize = true;

    // Create a simplified 1D mesh and attempt to re-create boundary conditions
    meshPtr.reset(new Foam::simplifiedMeshes::columnFvMesh(runTime));

    // Stopping after 1 iteration of the simplified mesh
    // Note: using saNoWriteNow will only trigger the function object execute
    // function and not the write function
    runTime.stopAt(Foam::Time::saNoWriteNow);

    if (args.optionFound("dry-run-write"))
    {
        // Stopping after 1 iteration of the simplified mesh
        // Note: using saWriteNow to trigger writing/execution of function
        // objects
        runTime.stopAt(Foam::Time::saWriteNow);
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
