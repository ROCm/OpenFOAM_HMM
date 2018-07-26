// check the final time directory for the following:

// 1. volume fields
HashTable<word> volumeFields;

// 2. the fields for each cloud:
HashTable<HashTable<word>> cloudFields;

if (timeDirs.size())
{
    const word& lastTimeName = timeDirs.last().name();
    const fileName cloudPrefix(regionPrefix/cloud::prefix);

    IOobjectList objs(mesh, lastTimeName);

    forAllConstIters(objs, fieldIter)
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
    // Now check for lagrangian/<cloudName>
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

        IOobjectList cloudObjs
        (
            mesh,
            lastTimeName,
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

                const word& fieldName = obj->name();
                const word& fieldType = obj->headerClassName();

                if (cloudFieldTypes.found(fieldType))
                {
                    // Field name/type - ignore types that we don't handle
                    fieldsPerCloud.insert(fieldName, fieldType);
                }
            }
        }
    }

    // Only retain a cloud that actually has fields
    cloudFields.filterValues
    (
        [](const HashTable<word>& v){ return v.size(); }
    );


    //
    // Verify that the variable is present for all times
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

        for (const fileName& file : contents)
        {
            missing.erase(file.name());
        }

        volumeFields.erase(missing);
    }
}


// ************************************************************************* //
