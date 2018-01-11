Foam::autoPtr<Foam::fvMesh> meshPtr(nullptr);

if (args.optionFound("check"))
{
    Info<< "Operating in 'check' mode: case will run for 1 time step.  "
        << "All checks assumed OK on a clean exit" << endl;

    Field<label>::allowConstructFromLargerSize = true;
    Field<scalar>::allowConstructFromLargerSize = true;
    Field<vector>::allowConstructFromLargerSize = true;
    Field<sphericalTensor>::allowConstructFromLargerSize = true;
    Field<symmTensor>::allowConstructFromLargerSize = true;
    Field<tensor>::allowConstructFromLargerSize = true;

    // Create a simplified 1D mesh and attempt to re-create boundary conditions
    meshPtr = dummyFvMesh::equivalent1DMesh(runTime);

    // Stopping after 1 iteration of the simplified mesh
    // Note: using saWriteNow to trigger writing/execution of function objects
    // - saNoWriteNow will only trigger the function object execute function
    runTime.stopAt(Foam::Time::saNoWriteNow);
//    runTime.stopAt(Foam::Time::saWriteNow);
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
