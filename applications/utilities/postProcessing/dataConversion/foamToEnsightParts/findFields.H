// check the final time directory for the following:

// 1. volume fields
HashTable<word> volumeFields;

// 2. the fields for each cloud:
HashTable<HashTable<word>> cloudFields;

if (timeDirs.size())
{
    const fileName& cloudPrefix = regionPrefix/cloud::prefix;
    const word& lastTimeName = timeDirs.last().name();

    IOobjectList objs(mesh, lastTimeName);

    forAllConstIter(IOobjectList, objs, fieldIter)
    {
        const IOobject& obj = *fieldIter();
        const word& fieldName = obj.name();
        const word& fieldType = obj.headerClassName();

        if (volFieldTypes.found(fieldType) && !fieldName.endsWith("_0"))
        {
            // ignore types that we don't handle, and ignore _0 fields
            volumeFields.insert(fieldName, fieldType);
        }
    }


    //
    // now check for lagrangian/<cloudName>
    //
    fileNameList cloudDirs;
    if (!noLagrangian)
    {
        cloudDirs = readDir
        (
            runTime.path()
          / lastTimeName
          / cloudPrefix,
            fileName::DIRECTORY
        );
    }

    forAll(cloudDirs, cloudI)
    {
        const word& cloudName = cloudDirs[cloudI];

        // Create a new hash table for each cloud
        cloudFields.insert(cloudName, HashTable<word>());

        // Identify the new cloud within the hash table
        HashTable<HashTable<word>>::iterator cloudIter =
            cloudFields.find(cloudName);

        IOobjectList objs
        (
            mesh,
            lastTimeName,
            cloudPrefix/cloudName
        );

        bool hasCoordinates = false;
        forAllConstIter(IOobjectList, objs, fieldIter)
        {
            const IOobject obj = *fieldIter();
            const word& fieldName = obj.name();
            const word& fieldType = obj.headerClassName();

            if (fieldName == "positions" || fieldName == "coordinates")
            {
                hasCoordinates = true;
            }
            else if (cloudFieldTypes.found(fieldType))
            {
                // simply ignore types that we don't handle
                cloudIter().insert(fieldName, fieldType);
            }
        }

        // drop this cloud if it has no positions or is otherwise empty
        if (!hasCoordinates || cloudIter().empty())
        {
            Info<< "removing cloud " << cloudName << endl;
            cloudFields.erase(cloudIter);
        }
    }

    //
    // verify that the variable is present for all times
    //
    for (label i=0; volumeFields.size() && i < timeDirs.size(); ++i)
    {
        const word& timeName = timeDirs[i].name();

        // Everything is potentially missing, unless we discover otherwise
        wordHashSet missing(volumeFields);

        // Avoid  -->>  IOobjectList objs(mesh, timeName);  <<--
        // Too much overhead when done so frequently.

        fileNameList contents = readDir
        (
            runTime.path()
          / timeName,
            fileName::FILE
        );

        forAll(contents, fileI)
        {
            missing.erase(contents[fileI].name());
        }

        volumeFields.erase(missing);
    }
}
