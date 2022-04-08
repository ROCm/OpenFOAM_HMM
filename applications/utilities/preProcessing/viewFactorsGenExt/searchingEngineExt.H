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

const triSurface localSurface = triangulate
(
    patches,
    includePatches,
    finalAgglom,
    triSurfaceToAgglom,
    globalNumbering,
    coarsePatches
);


distributedTriSurfaceMesh surfacesMesh
(
    IOobject
    (
        "wallSurface.stl",
        runTime.constant(),     // directory
        "triSurface",           // instance
        runTime,                // registry
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    localSurface,
    dict
);


triSurfaceToAgglom.resize(surfacesMesh.size());

// pbrt surface
std::vector<Point3f> vertices;

const pointField& pts = surfacesMesh.localPoints();
for (const auto& pt : pts)
{
    vertices.push_back(Point3f(pt[0], pt[1], pt[2]));
}

std::vector<int> indices;
for (const auto& tri : surfacesMesh.localFaces())
{
    indices.push_back(tri[0]);
    indices.push_back(tri[1]);
    indices.push_back(tri[2]);
}

std::vector<int> faceIndices;
faceIndices.reserve(surfacesMesh.localFaces().size());

for (label faceI = 0; faceI < surfacesMesh.localFaces().size(); faceI++)
{
    faceIndices.push_back(faceI);
}

Transform o2w;
Transform w2o;

std::vector<std::shared_ptr<Shape>> surfacesMesh_pbrt
(
    CreateTriangleMesh
    (
        &o2w,
        &w2o,
        false,                  // bool reverseOrientation
        surfacesMesh.size(),    // int nTriangles
        &indices[0],            // int *vertexIndices
        vertices.size(),
        &vertices[0],
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        nullptr,
        &faceIndices[0]
    )
);

std::vector<std::shared_ptr<Primitive>> prims;

prims.reserve(surfacesMesh_pbrt.size());
for (auto s : surfacesMesh_pbrt)
{
    prims.push_back
    (
        std::make_shared<GeometricPrimitive>
            (s, nullptr, nullptr, nullptr)
    );
}

std::shared_ptr<Primitive> accel;

ParamSet paramSet;

accel = CreateBVHAccelerator(std::move(prims), paramSet);
