Random rndGen(653213);

// Determine mesh bounding boxes:
List<treeBoundBox> meshBb
(
    1,
    treeBoundBox
    (
        boundBox(coarseMesh.points(), false)
    ).extend(rndGen, 1e-3)
);

// Dummy bounds dictionary
dictionary dict;
dict.add("bounds", meshBb);
dict.add
(
    "distributionType",
    distributedTriSurfaceMesh::distributionTypeNames_
    [
        distributedTriSurfaceMesh::FROZEN
    ]
);
dict.add("mergeDistance", SMALL);

labelList triSurfaceToAgglom(5*nFineFaces);

triSurface localSurface = triangulate
(
    patches,
    includePatches,
    finalAgglom,
    triSurfaceToAgglom,
    globalNumbering,
    coarsePatches
);


// CGAL surface

distributedTriSurfaceMesh surfacesMesh
(
    IOobject
    (
        "wallSurface.stl",
        runTime.constant(),     // directory
        "triSurface",           // instance
        runTime,                // registry
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    localSurface,
    dict
);

triSurfaceToAgglom.resize(surfacesMesh.size());

surfacesMesh.setField(triSurfaceToAgglom);

const pointField& pts = surfacesMesh.localPoints();

std::list<Triangle> triangles;

for (const auto& triLocal : surfacesMesh.localFaces())
{
    point p1l = pts[triLocal[0]];
    point p2l = pts[triLocal[1]];
    point p3l = pts[triLocal[2]];

    Point p1(p1l[0], p1l[1], p1l[2]);
    Point p2(p2l[0], p2l[1], p2l[2]);
    Point p3(p3l[0], p3l[1], p3l[2]);

    Triangle tri(p1, p2, p3);

    if (tri.is_degenerate())
    {
        std::cout << tri << std::endl;
    }
    triangles.push_back(tri);
}

// constructs AABB tree
Tree tree(triangles.begin(), triangles.end());
