// ignore special fields or fields that we don't handle
//
forAllIter(HashTable<word>, volumeFields, fieldIter)
{
    const word& fieldName = fieldIter.key();
    const word& fieldType = fieldIter();

    if (fieldName.size() > 2 &&  fieldName(fieldName.size() - 2, 2) != "_0")
    {
        volumeFields.erase(fieldIter);
    }
    else
    {
        // remove types that we don't handle:
        bool invalid = true;
        for (label typeI=0; invalid && volFieldTypes[typeI].size(); ++typeI)
        {
            if (fieldType == volFieldTypes[typeI])
            {
                invalid = false;
            }
        }

        if (invalid)
        {
            Info<< "ignoring " << fieldType << ": " << fieldName << endl;
            volumeFields.erase(fieldIter);
        }
    }
}

// verify that the variable is present for all times
//
for (label i=0; volumeFields.size() && i < timeDirs.size(); ++i)
{
    IOobjectList objs(mesh, timeDirs[i].name());

    forAllIter(HashTable<word>, volumeFields, fieldIter)
    {
        const word& fieldName = fieldIter.key();

        if (!objs.found(fieldName))
        {
            volumeFields.erase(fieldIter);
        }
    }
}



// ignore fields that we don't handle
//
forAllIter(HashTable<HashTable<word> >, cloudFields, cloudIter)
{
    const word& cloudName = cloudIter.key();

    forAllIter(HashTable<word>, cloudIter(), fieldIter)
    {
        const word& fieldName = fieldIter.key();
        const word& fieldType = fieldIter();

        // remove types that we don't handle:
        bool invalid = true;
        for (label typeI=0; invalid && sprayFieldTypes[typeI].size(); ++typeI)
        {
            if (fieldType == sprayFieldTypes[typeI])
            {
                invalid = false;
            }
        }

        if (invalid)
        {
            Info<< "ignoring " << fieldType << ": " << fieldName << endl;
            cloudIter().erase(fieldIter);
        }
    }

    if (!cloudIter().size())
    {
        Info<< "removing cloud " << cloudName<< endl;
        cloudFields.erase(cloudIter);
    }
}

// DEBUGGING
// Info<<"final fields (";
// forAllConstIter(HashTable<word>, volumeFields, fieldIter)
// {
//     Info<< " " << fieldIter.key();
// }
//
// Info<< " )\n";
//
// forAllConstIter(HashTable<HashTable<word> >, cloudFields, cloudIter)
// {
//     Info<<"final fields for lagrangian/" << cloudIter.key() << " (";
//     forAllConstIter(HashTable<word>, cloudIter(), fieldIter)
//     {
//         Info<< " " << fieldIter.key();
//     }
//     Info<< " )\n";
// }

