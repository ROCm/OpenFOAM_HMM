// check all time directories for the following:

// Any cloud names:
HashSet<word> allCloudDirs;

if (timeDirs.size() && !noLagrangian)
{
    const fileName& baseDir = mesh.time().path();
    const fileName cloudPrefix(regionPrefix/cloud::prefix);

    Info<< "Searching for lagrangian ... " << flush;

    forAll(timeDirs, timeI)
    {
        const word& timeName = timeDirs[timeI].name();

        // DO NOT USE -->>  runTime.setTime(timeDirs[timeI], timeI);  <<--
        // It incurs a large overhead when done so frequently.

        fileNameList cloudDirs = readDir
        (
            baseDir/timeName/cloudPrefix,
            fileName::DIRECTORY
        );

        for (const word& cloudName : cloudDirs)
        {
            IOobjectList cloudObjs
            (
                mesh,
                timeName,
                cloudPrefix/cloudName
            );

            // Clouds require "coordinates".
            // The "positions" are for v1706 and lower.
            if (cloudObjs.found("coordinates") || cloudObjs.found("positions"))
            {
                if (allCloudDirs.insert(cloudName))
                {
                    Info<< nl << "    At time: " << timeName
                        << " detected cloud directory : " << cloudName
                        << flush;
                }
            }
        }
    }

    if (allCloudDirs.empty())
    {
        Info<< "none detected." << endl;
    }

    if (Pstream::parRun())
    {
        Pstream::combineGather(allCloudDirs, HashSetPlusEqOp<word>());
        Pstream::combineScatter(allCloudDirs);
    }
}


// Sorted list of cloud names
const wordList cloudNames(allCloudDirs.sortedToc());

if (cloudNames.size())
{
    // Complete the echo information
    Info<< "(";
    for (const word& cloudName : cloudNames)
    {
        Info<< ' ' << cloudName;
    }
    Info<< " ) " << endl;
}


// ************************************************************************* //
