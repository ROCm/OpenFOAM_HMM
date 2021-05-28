/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Code chunk for post-processing conversion of cloud(s) to Ensight

\*---------------------------------------------------------------------------*/

// Cloud field data output
if (doLagrangian)
{
    // Lagrangian
    const auto& cloudFields = regionCloudFields[regioni];
    const auto& cloudNames = regionCloudNames[regioni];

    for (const word& cloudName : cloudNames)
    {
        const HashTable<word>& theseCloudFields = cloudFields[cloudName];

        fileNameList currentCloudDirs
        (
            readDir
            (
                runTime.timePath()/regionDir/cloud::prefix,
                fileName::DIRECTORY
            )
        );

        Info<< "Cloud " << cloudName << " (";

        const bool cloudExists =
            returnReduce(currentCloudDirs.found(cloudName), orOp<bool>());

        {
            autoPtr<ensightFile> os = ensCase.newCloud(cloudName);

            ensightOutput::writeCloudPositions
            (
                mesh,
                cloudName,
                cloudExists,
                os
            );

            Info<< " positions";
            if (!cloudExists)
            {
                Info<< "{0}"; // report empty field
            }
        }

        // Field order may differ on individual processors, so sort by name
        for (const word& fieldName : theseCloudFields.sortedToc())
        {
            const word& fieldType = theseCloudFields[fieldName];

            IOobject fieldObject
            (
                fieldName,
                mesh.time().timeName(),
                cloud::prefix/cloudName,
                mesh,
                IOobject::MUST_READ
            );

            bool fieldExists = cloudExists; // No field without positions
            if (cloudExists)
            {
                // Want MUST_READ (globally) and valid=false (locally),
                // but that combination does not work.
                // So check the header and sync globally

                const bool oldParRun = Pstream::parRun(false);
                fieldExists = fieldObject.typeHeaderOk<IOField<scalar>>(false);
                Pstream::parRun(oldParRun);  // Restore parallel state
                reduce(fieldExists, orOp<bool>());
            }

            bool wrote = false;
            if (fieldType == scalarIOField::typeName)
            {
                autoPtr<ensightFile> os =
                    ensCase.newCloudData<scalar>(cloudName, fieldName);

                wrote = ensightOutput::writeCloudField<scalar>
                (
                    fieldObject, fieldExists, os
                );
            }
            else if (fieldType == vectorIOField::typeName)
            {
                autoPtr<ensightFile> os =
                    ensCase.newCloudData<vector>(cloudName, fieldName);

                wrote = ensightOutput::writeCloudField<vector>
                (
                    fieldObject, fieldExists, os
                );
            }

            if (wrote)
            {
                Info<< ' ' << fieldName;
                if (!fieldExists)
                {
                    Info<< "{0}"; // report empty field
                }
            }
        }
        Info<< " )" << nl;
    }
}


// ************************************************************************* //
