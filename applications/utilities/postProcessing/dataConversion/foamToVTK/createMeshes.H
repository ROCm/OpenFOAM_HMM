PtrList<fvMesh> meshes(regionNames.size());
PtrList<fvMeshSubsetProxy> meshProxies(regionNames.size());
PtrList<vtk::vtuCells> vtuMappings(regionNames.size());

forAll(regionNames, regioni)
{
    const word& regionName = regionNames[regioni];

    Info<< "Create mesh";

    if (regionName != fvMesh::defaultRegion)
    {
        Info<< ' ' << regionName;
    }
    Info<< " for time = " << runTime.timeName() << nl << endl;

    meshes.set
    (
        regioni,
        new fvMesh
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

    // Mesh subsetting, or pass through
    meshProxies.set
    (
        regioni,
        new fvMeshSubsetProxy
        (
            meshes[regioni],
            cellSubsetType,
            cellSelectionName
        )
    );

    // VTU sizing and decomposition information
    vtuMappings.set
    (
        regioni,
        new vtk::vtuCells(writeOpts, decomposePoly)
    );
}


Info<< "VTK mesh topology: "
    << timer.cpuTimeIncrement() << " s, "
    << mem.update().size() << " kB" << endl;


// ************************************************************************* //
