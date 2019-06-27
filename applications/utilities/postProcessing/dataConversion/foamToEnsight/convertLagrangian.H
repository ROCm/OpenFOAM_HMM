/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Code chunk for post-processing conversion of cloud(s) to Ensight

\*---------------------------------------------------------------------------*/

// Cloud field data output
if (doLagrangian)
{
    forAll(cloudNames, cloudNo)
    {
        const word& cloudName = cloudNames[cloudNo];
        const HashTable<word>& theseCloudFields = cloudFields[cloudName];

        fileNameList currentCloudDirs
        (
            readDir
            (
                runTime.timePath()/regionPrefix/cloud::prefix,
                fileName::DIRECTORY
            )
        );

        Info<< "Write " << cloudName << " (";

        const bool cloudExists =
            returnReduce
            (
                currentCloudDirs.found(cloudName),
                orOp<bool>()
            );


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

        forAllConstIters(theseCloudFields, fieldIter)
        {
            const word& fieldName = fieldIter.key();
            const word& fieldType = fieldIter.val();

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

                fieldExists =
                    fieldObject.typeHeaderOk<IOField<scalar>>(false);

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
