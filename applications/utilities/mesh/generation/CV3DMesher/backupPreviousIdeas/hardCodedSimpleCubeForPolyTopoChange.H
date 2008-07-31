    IOobject io
    (
        Foam::polyMesh::defaultRegion,
        runTime.constant(),
        runTime,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    );

    polyMesh pMesh(io);

    IOobject io2
    (
        Foam::polyMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    );

    polyMesh pMesh2(io2, pointField(0), faceList(0), cellList(0));

    polyTopoChange meshCreator(pMesh.boundaryMesh().size());

    meshCreator.addPoint(point(-1,-1,-1), -1, -1, true);
    meshCreator.addPoint(point(1,-1,-1), -1, -1, true);
    meshCreator.addPoint(point(1,1,-1), -1, -1, true);
    meshCreator.addPoint(point(-1,1,-1), -1, -1, true);
    meshCreator.addPoint(point(-1,-1,1), -1, -1, true);
    meshCreator.addPoint(point(1,-1,1), -1, -1, true);
    meshCreator.addPoint(point(1,1,1), -1, -1, true);
    meshCreator.addPoint(point(-1,1,1), -1, -1, true);
    
    meshCreator.addCell(-1, -1, -1, -1, -1);

    labelList faceConstruct(4);

    faceConstruct[0] = 1;
    faceConstruct[1] = 2;
    faceConstruct[2] = 6;
    faceConstruct[3] = 5;

    meshCreator.addFace(face(faceConstruct), 0, -1, -1, -1, -1, false, 0, -1, false);

    faceConstruct[0] = 0;
    faceConstruct[1] = 4;
    faceConstruct[2] = 7;
    faceConstruct[3] = 3;

    meshCreator.addFace(face(faceConstruct), 0, -1, -1, -1, -1, false, 0, -1, false);

    faceConstruct[0] = 2;
    faceConstruct[1] = 3;
    faceConstruct[2] = 7;
    faceConstruct[3] = 6;

    meshCreator.addFace(face(faceConstruct), 0, -1, -1, -1, -1, false, 0, -1, false);

    faceConstruct[0] = 0;
    faceConstruct[1] = 1;
    faceConstruct[2] = 5;
    faceConstruct[3] = 4;

    meshCreator.addFace(face(faceConstruct), 0, -1, -1, -1, -1, false, 0, -1, false);

    faceConstruct[0] = 0;
    faceConstruct[1] = 3;
    faceConstruct[2] = 2;
    faceConstruct[3] = 1;

    meshCreator.addFace(face(faceConstruct), 0, -1, -1, -1, -1, false, 0, -1, false);

    faceConstruct[0] = 4;
    faceConstruct[1] = 5;
    faceConstruct[2] = 6;
    faceConstruct[3] = 7;

    meshCreator.addFace(face(faceConstruct), 0, -1, -1, -1, -1, false, 0, -1, false);

    Info<< meshCreator.points() << endl;
    Info<< meshCreator.faces() << endl;
    Info<< meshCreator.faceOwner() << endl;
    Info<< meshCreator.faceNeighbour() << endl;

//     calcDualMesh(meshCreator);

    autoPtr<mapPolyMesh> map = meshCreator.changeMesh(pMesh, false);

    pMesh.updateMesh(map);

    if (!pMesh.write())
    {
        FatalErrorIn("CV3D::writeMesh(const Time& runTime)")
            << "Failed writing polyMesh."
            << exit(FatalError);
    }