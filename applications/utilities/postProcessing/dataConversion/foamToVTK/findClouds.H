// check all time directories for the following:

// Any cloud names:
HashSet<fileName> allCloudDirs;

if (timeDirs.size() && !noLagrangian)
{
    const fileName& baseDir = mesh.time().path();
    const fileName& cloudPrefix = regionPrefix/cloud::prefix;

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

            // Clouds always require "positions"/"coordinates"
            if (cloudObjs.found("positions") || cloudObjs.found("coordinates"))
            {
                if (allCloudDirs.insert(cloudName))
                {
                    Info<< "At time: " << timeName
                        << " detected cloud directory : " << cloudName
                        << endl;
                }
            }
        }
    }

    if (allCloudDirs.empty())
    {
        Info<< "none detected." << endl;
    }
}

// sorted list of cloud names
const fileNameList cloudNames(allCloudDirs.sortedToc());

if (cloudNames.size())
{
    // complete the echo information
    Info<< "(";
    for (const word& cloudName : cloudNames)
    {
        Info<< ' ' << cloudName;
    }
    Info<< " ) " << endl;
}


// ************************************************************************* //
