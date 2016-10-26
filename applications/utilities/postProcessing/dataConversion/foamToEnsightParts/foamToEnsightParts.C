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

      - \par -index \<start\>
        Ignore the time index contained in the time file and use a
        simple indexing when creating the \c Ensight/data/######## files.

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

#include "ensightParts.H"
#include "ensightOutputFunctions.H"

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
        "define sub-directory name to use for Ensight data "
        "(default: \"Ensight\")"
    );
    argList::addOption
    (
        "width",
        "n",
        "width of Ensight data subdir"
    );

    // The volume field types that we handle
    wordHashSet volFieldTypes;
    volFieldTypes.insert(volScalarField::typeName);
    volFieldTypes.insert(volVectorField::typeName);
    volFieldTypes.insert(volSphericalTensorField::typeName);
    volFieldTypes.insert(volSymmTensorField::typeName);
    volFieldTypes.insert(volTensorField::typeName);

    // The lagrangian field types that we handle
    wordHashSet cloudFieldTypes;
    cloudFieldTypes.insert(scalarIOField::typeName);
    cloudFieldTypes.insert(vectorIOField::typeName);
    cloudFieldTypes.insert(tensorIOField::typeName);

    const char* geometryName = "geometry";

    #include "setRootCase.H"

    cpuTime timer;
    memInfo mem;
    Info<< "Initial memory "
        << mem.update().size() << " kB" << endl;

    #include "createTime.H"

    // Get times list
    instantList timeDirs = timeSelector::select0(runTime, args);

    // Default to binary output, unless otherwise specified
    const IOstream::streamFormat format =
    (
        args.optionFound("ascii")
      ? IOstream::ASCII
      : IOstream::BINARY
    );

    // Control for renumbering iterations
    label indexingNumber = 0;
    const bool optIndex = args.optionReadIfPresent("index", indexingNumber);
    const bool noLagrangian = args.optionFound("noLagrangian");

    // Always write the geometry, unless the -noMesh option is specified
    bool optNoMesh = args.optionFound("noMesh");

    // Adjust output width
    if (args.optionFound("width"))
    {
        ensightFile::subDirWidth(args.optionRead<label>("width"));
    }

    // Define sub-directory name to use for Ensight data
    fileName ensightDir = "Ensight";
    args.optionReadIfPresent("name", ensightDir);

    if (!ensightDir.isAbsolute())
    {
        ensightDir = args.rootPath()/args.globalCaseName()/ensightDir;
    }

    const fileName caseFileName = "Ensight.case";
    const fileName dataDir  = ensightDir/"data";
    const fileName dataMask = dataDir.name()/ensightFile::mask();

    // Ensight and Ensight/data directories must exist
    // do not remove old data - we might wish to convert new results
    // or a particular time interval
    if (isDir(ensightDir))
    {
        Info<<"Warning: re-using existing directory" << nl
            << "    " << ensightDir << endl;
    }

    // As per mkdir -p "Ensight/data"
    mkDir(ensightDir);
    mkDir(dataDir);

    #include "createNamedMesh.H"

    // Mesh instance (region0 gets filtered out)
    fileName regionPrefix;

    if (regionName != polyMesh::defaultRegion)
    {
        regionPrefix = regionName;
    }

    if (Pstream::master())
    {
        Info<< "Converting " << timeDirs.size() << " time steps" << endl;
    }

    // Construct the list of ensight parts for the entire mesh
    ensightParts partsList(mesh);

    // Write summary information
    {
        OFstream partsInfoFile(ensightDir/"partsInfo");

        partsInfoFile
            << "// summary of ensight parts" << nl << nl;
        partsList.writeSummary(partsInfoFile);
    }

    #include "checkHasMovingMesh.H"
    #include "findFields.H"

    if (hasMovingMesh && optNoMesh)
    {
        Info<< "mesh is moving: ignoring '-noMesh' option" << endl;
        optNoMesh = false;
    }


    // Map times used
    Map<scalar>  timeIndices;

    // TODO: Track the time indices used by the geometry
    DynamicList<label> geometryTimesUsed;

    // Track the time indices used by the volume fields
    DynamicList<label> fieldTimesUsed;

    // Track the time indices used by each cloud
    HashTable<DynamicList<label>> cloudTimesUsed;

    // Create a new DynamicList for each cloud
    forAllConstIter(HashTable<HashTable<word>>, cloudFields, cloudIter)
    {
        cloudTimesUsed.insert(cloudIter.key(), DynamicList<label>());
    }

    Info<< "Startup in "
        << timer.cpuTimeIncrement() << " s, "
        << mem.update().size() << " kB" << nl << endl;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        #include "getTimeIndex.H"

        // Remember the time index for the volume fields
        fieldTimesUsed.append(timeIndex);

        // The data/ITER subdirectory must exist
        // Note that data/ITER is indeed a valid ensight::FileName
        const fileName subDir = ensightFile::subDir(timeIndex);
        mkDir(dataDir/subDir);

        // Place a timestamp in the directory for future reference
        {
            OFstream timeStamp(dataDir/subDir/"time");
            timeStamp
                << "#   timestep time" << nl
                << subDir.c_str() << " " << runTime.timeName() << nl;
        }

        #include "moveMesh.H"

        if (timeI == 0 || mesh.moving())
        {
            if (mesh.moving())
            {
                partsList.recalculate(mesh);
            }

            if (!optNoMesh)
            {
                if (hasMovingMesh)
                {
                    // Remember the time index for the geometry
                    geometryTimesUsed.append(timeIndex);
                }

                ensightGeoFile geoFile
                (
                    (hasMovingMesh ? dataDir/subDir : ensightDir),
                    geometryName,
                    format
                );
                partsList.writeGeometry(geoFile);
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

            if (fieldType == volScalarField::typeName)
            {
                ensightVolField<scalar>
                (
                    partsList,
                    fieldObject,
                    mesh,
                    dataDir,
                    subDir,
                    format
                );

            }
            else if (fieldType == volVectorField::typeName)
            {
                ensightVolField<vector>
                (
                    partsList,
                    fieldObject,
                    mesh,
                    dataDir,
                    subDir,
                    format
                );

            }
            else if (fieldType == volSphericalTensorField::typeName)
            {
                ensightVolField<sphericalTensor>
                (
                    partsList,
                    fieldObject,
                    mesh,
                    dataDir,
                    subDir,
                    format
                );

            }
            else if (fieldType == volSymmTensorField::typeName)
            {
                ensightVolField<symmTensor>
                (
                    partsList,
                    fieldObject,
                    mesh,
                    dataDir,
                    subDir,
                    format
                );
            }
            else if (fieldType == volTensorField::typeName)
            {
                ensightVolField<tensor>
                (
                    partsList,
                    fieldObject,
                    mesh,
                    dataDir,
                    subDir,
                    format
                );
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

            // Check that the positions field is present for this time
            if (!cloudObjs.found("positions"))
            {
                continue;
            }

            Info<< "Write " << cloudName << " ( positions" << flush;

            ensightParticlePositions
            (
                mesh,
                dataDir,
                subDir,
                cloudName,
                format
            );

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

                if (fieldType == scalarIOField::typeName)
                {
                    ensightLagrangianField<scalar>
                    (
                        *fieldObject,
                        dataDir,
                        subDir,
                        cloudName,
                        format
                    );

                }
                else if (fieldType == vectorIOField::typeName)
                {
                    ensightLagrangianField<vector>
                    (
                        *fieldObject,
                        dataDir,
                        subDir,
                        cloudName,
                        format
                    );

                }
                else if (fieldType == tensorIOField::typeName)
                {
                    ensightLagrangianField<tensor>
                    (
                        *fieldObject,
                        dataDir,
                        subDir,
                        cloudName,
                        format
                    );

                }
            }

            Info<< " )" << endl;

            // Remember the time index
            cloudTimesUsed[cloudName].append(timeIndex);
        }

        Info<< "Wrote in "
            << timer.cpuTimeIncrement() << " s, "
            << mem.update().size() << " kB" << endl;
    }

    #include "ensightOutputCase.H"

    Info<< "\nEnd: "
        << timer.elapsedCpuTime() << " s, "
        << mem.update().peak() << " kB (peak)\n" << endl;

    return 0;
}


// ************************************************************************* //
