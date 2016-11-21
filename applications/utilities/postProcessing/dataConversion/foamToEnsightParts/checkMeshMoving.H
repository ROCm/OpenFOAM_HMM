// check for "points" in all of the result directories
// - could restrict to the selected times

bool meshMoving = false;

if (timeDirs.size() > 1 && Pstream::master())
{
    // We already loaded a mesh (usually from constant).
    // See if any other "polyMesh/points" files exist too.

    Info<< "Search for moving mesh ... " << flush;
    forAll(timeDirs, timeI)
    {
        meshMoving =
        (
            IOobject
            (
                "points",
                timeDirs[timeI].name(),
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
        Info<< "found." << nl
            << "    Writing meshes for every timestep." << endl;
    }
    else
    {
        Info<< "none detected." << endl;
    }
}

reduce(meshMoving, orOp<bool>());
