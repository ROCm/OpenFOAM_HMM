// check for "points" in any of the result directories

bool meshMoving = false;
if (timeDirs.size() > 1)
{
    // We already loaded a mesh (usually from constant). See if any other
    // points files
    forAll(timeDirs, timeI)
    {
        if (timeDirs[timeI].name() != mesh.pointsInstance())
        {
            meshMoving = IOobject
            (
                "points",
                timeDirs[timeI].name(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::NO_READ
            ).typeHeaderOk<pointIOField>(true);

            if (meshMoving)
            {
                break;
            }
        }
    }
}
