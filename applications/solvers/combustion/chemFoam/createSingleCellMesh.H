Info<< "Constructing single cell mesh" << nl << endl;

autoPtr<fvMesh> meshPtr = dummyFvMesh::singleCellMesh(runTime);
fvMesh& mesh = meshPtr();
