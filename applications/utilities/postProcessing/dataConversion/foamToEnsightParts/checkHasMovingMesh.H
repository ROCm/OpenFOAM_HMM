// check for "points" in all of the result directories
// - could restrict to the selected times

bool hasMovingMesh = false;

if (timeDirs.size() > 1 && Pstream::master())
{
    // We already loaded a mesh (usually from constant).
    // See if any other "polyMesh/points" files exist too.

    const fileName& baseDir = mesh.time().path();

    Info<< "Search for moving mesh ... " << flush;
    forAll(timeDirs, timeI)
    {
        hasMovingMesh =
        (
            isDir(baseDir/timeDirs[timeI].name()/polyMesh::meshSubDir)
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

        if (hasMovingMesh)
        {
            break;
        }
    }

    if (hasMovingMesh)
    {
        Info<< "found." << nl
            << "    Writing meshes for every timestep." << endl;
    }
    else
    {
        Info<< "none detected." << endl;
    }
}

reduce(hasMovingMesh, orOp<bool>());
