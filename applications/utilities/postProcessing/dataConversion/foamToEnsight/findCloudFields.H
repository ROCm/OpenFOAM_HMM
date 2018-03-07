// check all time directories for the following:

// The fields for each cloud:
HashTable<HashTable<word>> cloudFields;

// Identify if lagrangian data exist at any time step.
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

        forAll(cloudDirs, cloudI)
        {
            const word& cloudName = cloudDirs[cloudI];

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
                // Save the cloud fields on a per cloud basis
                auto& fieldsPerCloud = cloudFields(cloudName);

                forAllConstIters(cloudObjs, fieldIter)
                {
                    const IOobject* obj = fieldIter();

                    // Field name/type
                    fieldsPerCloud.insert(obj->name(), obj->headerClassName());
                }
            }
        }
    }

    // Prune out geometry again since it gets treated specially
    forAllIters(cloudFields, cloudIter)
    {
        cloudIter().erase("coordinates");
        cloudIter().erase("positions");
    }

    if (Pstream::parRun())
    {
        Pstream::mapCombineGather(cloudFields, HashTableOps::plusEqOp<word>());
        Pstream::mapCombineScatter(cloudFields);
    }

    if (cloudFields.empty())
    {
        Info<< "none detected." << endl;
    }
}


// Sorted list of cloud names
const wordList cloudNames(cloudFields.sortedToc());

if (cloudNames.size())
{
    // Complete the echo information - as flatOutput
    cloudNames.writeList(Info) << endl;
}


// ************************************************************************* //
