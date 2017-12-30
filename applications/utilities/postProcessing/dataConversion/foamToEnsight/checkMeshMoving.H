// Check for "points" in any of the result directories

bool meshMoving = false;

if (timeDirs.size() > 1 && Pstream::master())
{
    // We already loaded a mesh (usually from constant).
    // See if any other "polyMesh/points" files exist too.

    Info<< "Search for moving mesh ... " << flush;
    forAll(timeDirs, timeI)
    {
        const word& timeName = timeDirs[timeI].name();

        meshMoving =
        (
            timeName != mesh.pointsInstance()
         && IOobject
            (
                "points",
                timeName,
                polyMesh::meshSubDir,
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false  // no register
            ).typeHeaderOk<pointIOField>(true, false)
        );

        if (meshMoving)
        {
            break;
        }
    }

    if (meshMoving)
    {
        Info<< "found. Writing meshes for every timestep." << endl;
    }
    else
    {
        Info<< "none detected." << endl;
    }
}

reduce(meshMoving, orOp<bool>());
