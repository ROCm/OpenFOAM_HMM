// check all time directories for the following:

// The fields for each cloud:
HashTable<HashTable<word>> cloudFields;

// Identify if lagrangian data exist at any time step.
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

        forAll(cloudDirs, cloudI)
        {
            const word& cloudName = cloudDirs[cloudI];

            IOobjectList cloudObjs
            (
                mesh,
                timeName,
                cloudPrefix/cloudName
            );

            // clouds always require "positions"
            if (cloudObjs.found("positions"))
            {
                HashTable<HashTable<word>>::iterator cloudIter =
                    cloudFields.find(cloudName);

                if (cloudIter == cloudFields.end())
                {
                    // A newly discovered cloud
                    cloudFields.insert(cloudName, HashTable<word>());
                    cloudIter = cloudFields.find(cloudName);
                }

                forAllConstIter(IOobjectList, cloudObjs, fieldIter)
                {
                    const IOobject& obj = *fieldIter();

                    // Add field and field type
                    cloudIter().insert
                    (
                        obj.name(),
                        obj.headerClassName()
                    );
                }
            }
        }
    }

    // prune out "positions" again since it gets treated specially
    forAllIter(HashTable<HashTable<word>>, cloudFields, cloudIter)
    {
        cloudIter().erase("positions");
    }

    if (cloudFields.empty())
    {
        Info<< "none detected." << endl;
    }
}


// sorted list of cloud names
const wordList cloudNames(cloudFields.sortedToc());

if (cloudNames.size())
{
    // complete the echo information
    Info<< "(";
    forAll(cloudNames, cloudNo)
    {
        Info<< ' ' << cloudNames[cloudNo];
    }
    Info<< " ) " << endl;
}

// ************************************************************************* //
