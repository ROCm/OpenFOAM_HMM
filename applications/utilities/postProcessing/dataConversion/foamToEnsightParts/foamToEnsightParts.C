/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

Application
    foamToEnsightParts

Group
    grpPostProcessingUtilities

Description
    Translates OpenFOAM data to Ensight format.
    An Ensight part is created for each cellZone and patch.

Usage
    \b foamToEnsightParts [OPTION]

    Options:
      - \par -ascii
        Write Ensight data in ASCII format instead of "C Binary"

      - \par -name \<subdir\>
        Define sub-directory name to use for Ensight data (default: "Ensight")

      - \par -noZero
        Exclude the often incomplete initial conditions.

      - \par -index \<start\>
        Ignore the time index contained in the time file and use a
        simple indexing when creating the \c Ensight/data/######## files.

      - \par -noLagrangian
        Suppress writing lagrangian positions and fields.

      - \par -noMesh
        Suppress writing the geometry. Can be useful for converting partial
        results for a static geometry.

      - \par -width \<n\>
        Width of Ensight data subdir

Note
    - no parallel data.
    - writes to \a Ensight directory to avoid collisions with foamToEnsight.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"

#include "volFields.H"
#include "OFstream.H"
#include "IOmanip.H"
#include "IOobjectList.H"
#include "scalarIOField.H"
#include "tensorIOField.H"

// file-format/conversion
#include "ensightCase.H"
#include "ensightGeoFile.H"
#include "ensightParts.H"
#include "ensightSerialOutput.H"

// local files
#include "ensightOutputSerialCloud.H"

#include "memInfo.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    // Enable -constant
    // Probably don't need -withZero though, since the fields are vetted
    // afterwards anyhow
    timeSelector::addOptions(true, false);
    argList::noParallel();
    argList::addBoolOption
    (
        "ascii",
        "write in ASCII format instead of 'C Binary'"
    );
    argList::addOption
    (
        "index",
        "start",
        "ignore the time index contained in the uniform/time file "
        "and use simple indexing when creating the files"
    );
    argList::addBoolOption
    (
        "noLagrangian",
        "suppress writing lagrangian positions and fields"
    );
    argList::addBoolOption
    (
        "noMesh",
        "suppress writing the geometry. "
        "Can be useful for converting partial results for a static geometry"
    );
    argList::addOption
    (
        "name",
        "subdir",
        "sub-directory name for ensight output (default: 'Ensight')"
    );
    argList::addOption
    (
        "width",
        "n",
        "width of Ensight data subdir"
    );

    // The volume field types that we handle
    const wordHashSet volFieldTypes
    {
        volScalarField::typeName,
        volVectorField::typeName,
        volSphericalTensorField::typeName,
        volSymmTensorField::typeName,
        volTensorField::typeName
    };

    // The lagrangian field types that we handle
    const wordHashSet cloudFieldTypes
    {
        scalarIOField::typeName,
        vectorIOField::typeName,
        tensorIOField::typeName
    };

    #include "setRootCase.H"

    // Default to binary output, unless otherwise specified
    const IOstream::streamFormat format =
    (
        args.optionFound("ascii")
      ? IOstream::ASCII
      : IOstream::BINARY
    );

    cpuTime timer;
    memInfo mem;
    Info<< "Initial memory " << mem.update().size() << " kB" << endl;

    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createNamedMesh.H"

    fileName regionPrefix; // Mesh instance (region0 gets filtered out)
    if (regionName != polyMesh::defaultRegion)
    {
        regionPrefix = regionName;
    }

    //
    // general (case) output options
    //
    ensightCase::options caseOpts(format);

    caseOpts.width(args.optionLookupOrDefault<label>("width", 8));
    caseOpts.overwrite(false); // leave existing output directory

    // Can also have separate directory for lagrangian
    // caseOpts.separateCloud(true);

    // Define sub-directory name to use for EnSight data.
    // The path to the ensight directory is at case level only
    // - For parallel cases, data only written from master
    fileName ensightDir = args.optionLookupOrDefault<word>("name", "Ensight");
    if (!ensightDir.isAbsolute())
    {
        ensightDir = args.rootPath()/args.globalCaseName()/ensightDir;
    }

    //
    // Open new ensight case file, initialize header etc.
    //
    ensightCase ensCase
    (
        ensightDir,
        "Ensight",  // args.globalCaseName(),
        caseOpts
    );


    //
    // Miscellaneous output configuration
    //

    // Control for renumbering iterations
    label indexingNumber = 0;
    const bool optIndex = args.optionReadIfPresent("index", indexingNumber);
    const bool noLagrangian = args.optionFound("noLagrangian");

    // Always write the geometry, unless the -noMesh option is specified
    bool optNoMesh = args.optionFound("noMesh");


    // Construct the list of ensight parts for the entire mesh
    ensightParts partsList(mesh);

    // Write summary information
    if (Pstream::master())
    {
        Info<< "Converting " << timeDirs.size() << " time steps" << endl;

        OFstream info(ensCase.path()/"partsInfo");

        info
            << "// summary of ensight parts" << nl << nl;
        partsList.writeSummary(info);
    }

    #include "checkMeshMoving.H"
    #include "findFields.H"

    if (meshMoving && optNoMesh)
    {
        Info<< "mesh is moving: ignoring '-noMesh' option" << endl;
        optNoMesh = false;
    }


    Info<< "Startup in "
        << timer.cpuTimeIncrement() << " s, "
        << mem.update().size() << " kB" << nl << endl;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        #include "getTimeIndex.H"
        #include "moveMesh.H"

        ensCase.setTime(timeDirs[timeI], timeIndex);

        if (timeI == 0 || mesh.moving())
        {
            if (mesh.moving())
            {
                partsList.recalculate(mesh);
            }

            if (!optNoMesh)
            {
                autoPtr<ensightGeoFile> os = ensCase.newGeometry(meshMoving);
                partsList.write(os.rawRef());
            }
        }

        Info<< "Write volume field (" << flush;

        forAllConstIter(HashTable<word>, volumeFields, fieldIter)
        {
            const word& fieldName = fieldIter.key();
            const word& fieldType = fieldIter();

            IOobject fieldObject
            (
                fieldName,
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            );

            bool wrote = false;
            if (fieldType == volScalarField::typeName)
            {
                autoPtr<ensightFile> os = ensCase.newData<scalar>
                (
                    fieldName
                );

                volScalarField vf(fieldObject, mesh);
                wrote = ensightSerialOutput::writeField<scalar>
                (
                    vf, partsList, os
                );
            }
            else if (fieldType == volVectorField::typeName)
            {
                autoPtr<ensightFile> os = ensCase.newData<vector>
                (
                    fieldName
                );

                volVectorField vf(fieldObject, mesh);
                wrote = ensightSerialOutput::writeField<vector>
                (
                    vf, partsList, os
                );
            }
            else if (fieldType == volSphericalTensorField::typeName)
            {
                autoPtr<ensightFile> os = ensCase.newData<sphericalTensor>
                (
                    fieldName
                );

                volSphericalTensorField vf(fieldObject, mesh);
                wrote = ensightSerialOutput::writeField<sphericalTensor>
                (
                    vf, partsList, os
                );
            }
            else if (fieldType == volSymmTensorField::typeName)
            {
                autoPtr<ensightFile> os = ensCase.newData<symmTensor>
                (
                    fieldName
                );

                volSymmTensorField vf(fieldObject, mesh);
                wrote = ensightSerialOutput::writeField<symmTensor>
                (
                    vf, partsList, os
                );
            }
            else if (fieldType == volTensorField::typeName)
            {
                autoPtr<ensightFile> os = ensCase.newData<tensor>
                (
                    fieldName
                );

                volTensorField vf(fieldObject, mesh);
                wrote = ensightSerialOutput::writeField<tensor>
                (
                    vf, partsList, os
                );
            }

            if (wrote)
            {
                Info<< " " << fieldObject.name() << flush;
            }
        }
        Info<< " )" << endl;

        // Check for clouds
        forAllConstIter(HashTable<HashTable<word>>, cloudFields, cloudIter)
        {
            const word& cloudName = cloudIter.key();
            const fileName& cloudPrefix = regionPrefix/cloud::prefix;

            if (!isDir(runTime.timePath()/cloudPrefix/cloudName))
            {
                continue;
            }

            IOobjectList cloudObjs
            (
                mesh,
                runTime.timeName(),
                cloudPrefix/cloudName
            );

            // Check that the positions/coordinates field is present for this
            // time
            if
            (
                !cloudObjs.found("positions")
             || !cloudObjs.found("coordinates")
            )
            {
                continue;
            }

            Info<< "Write " << cloudName << " (" << flush;

            ensightSerialCloud::writePositions
            (
                mesh,
                cloudName,
                ensCase.newCloud(cloudName)
            );
            Info<< " positions";


            forAllConstIter(HashTable<word>, cloudIter(), fieldIter)
            {
                const word& fieldName = fieldIter.key();
                const word& fieldType = fieldIter();

                IOobject *fieldObject = cloudObjs.lookup(fieldName);

                if (!fieldObject)
                {
                    Info<< "missing "
                        << runTime.timeName()/cloudPrefix/cloudName
                        / fieldName
                        << endl;
                    continue;
                }

                bool wrote = false;
                if (fieldType == scalarIOField::typeName)
                {
                    wrote = ensightSerialCloud::writeCloudField<scalar>
                    (
                        *fieldObject,
                        ensCase.newCloudData<scalar>(cloudName, fieldName)
                    );
                }
                else if (fieldType == vectorIOField::typeName)
                {
                    wrote = ensightSerialCloud::writeCloudField<vector>
                    (
                        *fieldObject,
                        ensCase.newCloudData<vector>(cloudName, fieldName)
                    );
                }
                else if (fieldType == tensorIOField::typeName)
                {
                    wrote = ensightSerialCloud::writeCloudField<tensor>
                    (
                        *fieldObject,
                        ensCase.newCloudData<tensor>(cloudName, fieldName)
                    );
                }

                if (wrote)
                {
                    Info<< " " << fieldObject->name();
                }
            }

            Info<< " )" << endl;
        }

        Info<< "Wrote in "
            << timer.cpuTimeIncrement() << " s, "
            << mem.update().size() << " kB" << endl;
    }

    ensCase.write();

    Info<< "\nEnd: "
        << timer.elapsedCpuTime() << " s, "
        << mem.update().peak() << " kB (peak)\n" << endl;

    return 0;
}


// ************************************************************************* //
