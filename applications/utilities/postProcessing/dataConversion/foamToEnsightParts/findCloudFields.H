// check all time directories for the following:

// The fields for each cloud:
HashTable<HashTable<word>> cloudFields;

// Identify if lagrangian data exist at any time step.
if (timeDirs.size() && doLagrangian)
{
    const fileName& baseDir = mesh.time().path();
    const fileName cloudPrefix(regionPrefix/cloud::prefix);

    Info<< "Searching for lagrangian ... " << flush;

    for (const instant& inst : timeDirs)
    {
        const word& timeName = inst.name();

        // DO NOT USE -->>  runTime.setTime(timeDirs[timeI], timeI);  <<--
        // It incurs a large overhead when done so frequently.

        fileNameList cloudDirs
        (
            readDir
            (
                baseDir/timeName/cloudPrefix,
                fileName::DIRECTORY
            )
        );

        for (fileName& cloudDir : cloudDirs)
        {
            const word cloudName(std::move(cloudDir));

            IOobjectList cloudObjs
            (
                mesh,
                timeName,
                cloudPrefix/cloudName
            );

            // Clouds require "coordinates".
            // The "positions" are for v1706 and lower.
            // - detect and remove since these are treated specially

            bool isCloud = false;
            if (cloudObjs.erase("coordinates"))
            {
                isCloud = true;
            }
            if (cloudObjs.erase("positions"))
            {
                isCloud = true;
            }

            if (isCloud)
            {
                // Save the cloud fields on a per cloud basis
                auto& fieldsPerCloud = cloudFields(cloudName);

                forAllConstIters(cloudObjs, fieldIter)
                {
                    const IOobject* io = *fieldIter;

                    // Field name/type
                    fieldsPerCloud.insert(io->name(), io->headerClassName());
                }
            }
        }
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
