// check the final time directory for

// 1. volume fields
HashTable<word> volumeFields;

// 2. the fields for each cloud:
HashTable<HashTable<word> > cloudFields;

if (timeDirs.size() > 1)
{
    IOobjectList objs(mesh, timeDirs[timeDirs.size()-1].name());

    forAllConstIter(IOobjectList, objs, fieldIter)
    {
        const IOobject& obj = *fieldIter();

        if
        (
            obj.headerClassName() == volScalarField::typeName
         || obj.headerClassName() == volVectorField::typeName
        )
        {
            // Add field and field type
            volumeFields.insert
            (
                obj.name(),
                obj.headerClassName()
            );
        }
    }

    // now check for lagrangian/<cloudName>

    fileNameList cloudDirs = readDir
    (
        runTime.path()
      / timeDirs[timeDirs.size() - 1].name()
      / regionPrefix
      / "lagrangian",
        fileName::DIRECTORY
    );

    forAll(cloudDirs, cloudI)
    {
        const word& cloudName = cloudDirs[cloudI];

        // Create a new hash table for each cloud
        cloudFields.insert(cloudName, HashTable<word>());

        // Identify the new cloud in the hash table
        HashTable<HashTable<word> >::iterator cloudIter =
            cloudFields.find(cloudName);

        IOobjectList cloudObjs
        (
            mesh,
            timeDirs[timeDirs.size() - 1].name(),
            "lagrangian"/cloudName
        );

        bool hasPositions = false;
        forAllConstIter(IOobjectList, cloudObjs, fieldIter)
        {
            const IOobject obj = *fieldIter();

            if (obj.name() == "positions")
            {
                hasPositions = true;
            }
            else
            {
                // Add field and field type
                cloudIter().insert
                (
                    obj.name(),
                    obj.headerClassName()
                );
            }
        }

        // drop this cloud if it has no positions
        if (!hasPositions)
        {
            cloudFields.erase(cloudIter);
        }
    }
}

