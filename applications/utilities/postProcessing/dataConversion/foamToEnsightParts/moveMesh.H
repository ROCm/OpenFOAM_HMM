{
    IOobject io
    (
        "points",
        runTime.timeName(),
        polyMesh::meshSubDir,
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false  // no register
    );

    if (io.typeHeaderOk<pointIOField>(true, false))
    {
        // Read new points
        io.readOpt() = IOobject::MUST_READ;
        mesh.movePoints(pointIOField(io));
    }
}
