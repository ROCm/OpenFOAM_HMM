// check for "points" in any of the result directories

bool meshMoving = false;

if (timeDirs.size() > 1 && Pstream::master())
{
    // We already loaded a mesh (usually from constant).
    // See if any other "polyMesh/points" files exist too.

    const fileName& baseDir = mesh.time().path();

    Info<< "Search for moving mesh ... " << flush;
    forAll(timeDirs, timeI)
    {
        meshMoving =
        (
            timeDirs[timeI].name() != mesh.pointsInstance()
         && isDir(baseDir/timeDirs[timeI].name()/polyMesh::meshSubDir)
         && IOobject
            (
                "points",
                timeDirs[timeI].name(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false  // no register
            ).typeHeaderOk<pointIOField>(true)
        );

        if (meshMoving)
        {
            break;
        }
    }

    if (meshMoving)
    {
        Info<< "found." << nl
            << "    Writing meshes for every timestep." << endl;
    }
    else
    {
        Info<< "none detected." << endl;
    }
}

reduce(meshMoving, orOp<bool>());
