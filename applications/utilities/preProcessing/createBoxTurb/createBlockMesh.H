const cellModel& hex = cellModel::ref(cellModel::HEX);

cellShapeList cellShapes;
faceListList boundary;
pointField points;
{
    Info<< "Creating block" << endl;

    block b
    (
        cellShape(hex, identity(8)),
        pointField
        (
            {
                point(0, 0, 0),
                point(L.x(), 0, 0),
                point(L.x(), L.y(), 0),
                point(0, L.y(), 0),
                point(0, 0, L.z()),
                point(L.x(), 0, L.z()),
                point(L.x(), L.y(), L.z()),
                point(0, L.y(), L.z())
            }
        ),
        blockEdgeList(),
        blockFaceList(),
        N
    );

    Info<< "Creating cells" << endl;
    cellShapes = b.shapes();

    Info<< "Creating boundary faces" << endl;

    boundary.setSize(b.boundaryPatches().size());
    forAll(boundary, patchi)
    {
        faceList faces(b.boundaryPatches()[patchi].size());
        forAll(faces, facei)
        {
            faces[facei] = face(b.boundaryPatches()[patchi][facei]);
        }
        boundary[patchi].transfer(faces);
    }

    points.transfer(const_cast<pointField&>(b.points()));
}

Info<< "Creating patch dictionaries" << endl;
wordList patchNames(boundary.size());
forAll(patchNames, patchi)
{
    patchNames[patchi] = polyPatch::defaultName(patchi);
}

PtrList<dictionary> boundaryDicts(boundary.size());
forAll(boundaryDicts, patchi)
{
    boundaryDicts.set(patchi, new dictionary());
    dictionary& patchDict = boundaryDicts[patchi];
    word nbrPatchName;
    if (patchi % 2 == 0)
    {
        nbrPatchName = polyPatch::defaultName(patchi + 1);
    }
    else
    {
        nbrPatchName = polyPatch::defaultName(patchi - 1);
    }

    patchDict.add("type", cyclicPolyPatch::typeName);
    patchDict.add("neighbourPatch", nbrPatchName);
}

Info<< "Creating polyMesh" << endl;
polyMesh mesh
(
    IOobject
    (
        polyMesh::defaultRegion,
        runTime.constant(),
        runTime,
        IOobject::NO_READ
    ),
    std::move(points),
    cellShapes,
    boundary,
    patchNames,
    boundaryDicts,
    "defaultFaces",
    cyclicPolyPatch::typeName,
    false
);

Info<< "Writing polyMesh" << endl;
mesh.write();
