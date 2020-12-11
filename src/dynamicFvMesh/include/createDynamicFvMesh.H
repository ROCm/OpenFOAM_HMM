Info<< "Create mesh for time = "
    << runTime.timeName() << nl << endl;

autoPtr<dynamicFvMesh> meshPtr(dynamicFvMesh::New(args, runTime));//, false));

dynamicFvMesh& mesh = meshPtr();

// initialise all (lower levels and current)
mesh.init(true);
