/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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
    foamToEnsight

Group
    grpPostProcessingUtilitie

Description
    Translates OpenFOAM data to EnSight format.

    An Ensight part is created for the internalMesh and for each patch.

Usage
    \b foamToEnsight [OPTION]

    Options:
      - \par -ascii
        Write Ensight data in ASCII format instead of "C Binary"

      - \par -noZero
        Exclude the often incomplete initial conditions.

      - \par -noLagrangian
        Suppress writing lagrangian positions and fields.

      - \par -noPatches
        Suppress writing any patches.

      - \par -patches patchList
        Specify particular patches to write.
        Specifying an empty list suppresses writing the internalMesh.

      - \par -faceZones zoneList
        Specify faceZones to write, with wildcards

      - \par -cellZone zoneName
        Specify single cellZone to write (not lagrangian)

      - \par -width \<n\>
        Width of EnSight data subdir (default: 8)

Note
    Writes to \a EnSight directory to avoid collisions with
    foamToEnsightParts

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "IOobjectList.H"
#include "IOmanip.H"
#include "OFstream.H"

#include "fvc.H"
#include "volFields.H"
#include "hashedWordList.H"

#include "labelIOField.H"
#include "scalarIOField.H"
#include "tensorIOField.H"

// file-format/conversion
#include "ensightCase.H"
#include "ensightGeoFile.H"
#include "ensightMesh.H"
#include "ensightOutput.H"

// local files
#include "meshSubsetHelper.H"
#include "ensightOutputCloud.H"

#include "memInfo.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"

    argList::addBoolOption
    (
        "ascii",
        "write in ASCII format instead of 'C Binary'"
    );
    argList::addBoolOption
    (
        "nodeValues",
        "write values in nodes"
    );
    argList::addBoolOption
    (
        "noLagrangian",
        "suppress writing lagrangian positions and fields"
    );
    argList::addBoolOption
    (
        "noPatches",
        "suppress writing any patches"
    );
    argList::addOption
    (
        "patches",
        "wordReList",
        "specify particular patches to write - eg '(outlet \"inlet.*\")'. "
        "An empty list suppresses writing the internalMesh."
    );
    argList::addOption
    (
        "faceZones",
        "wordReList",
        "specify faceZones to write - eg '( slice \"mfp-.*\" )'."
    );
    argList::addOption
    (
        "fields",
        "wordReList",
        "specify fields to export (all by default) - eg '( \"U.*\" )'."
    );
    argList::addOption
    (
        "cellZone",
        "word",
        "specify cellZone to write"
    );
    argList::addOption
    (
        "name",
        "subdir",
        "sub-directory name for ensight output (default: 'EnSight')"
    );
    argList::addOption
    (
        "width",
        "n",
        "width of ensight data subdir"
    );

    // The volume field types that we handle
    const hashedWordList volFieldTypes
    {
        volScalarField::typeName,
        volVectorField::typeName,
        volSphericalTensorField::typeName,
        volSymmTensorField::typeName,
        volTensorField::typeName,

        volScalarField::Internal::typeName,
        volVectorField::Internal::typeName,
        volSphericalTensorField::Internal::typeName,
        volSymmTensorField::Internal::typeName,
        volTensorField::Internal::typeName
    };

    #include "setRootCase.H"

    // Default to binary output, unless otherwise specified
    const IOstream::streamFormat format =
    (
        args.optionFound("ascii")
      ? IOstream::ASCII
      : IOstream::BINARY
    );

    const bool nodeValues = args.optionFound("nodeValues");

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
    // General (case) output options
    //
    ensightCase::options caseOpts(format);

    caseOpts.nodeValues(args.optionFound("nodeValues"));
    caseOpts.width(args.optionLookupOrDefault<label>("width", 8));
    caseOpts.overwrite(true); // remove existing output directory

    // Can also have separate directory for lagrangian
    // caseOpts.separateCloud(true);


    // Define sub-directory name to use for EnSight data.
    // The path to the ensight directory is at case level only
    // - For parallel cases, data only written from master
    fileName ensightDir = args.optionLookupOrDefault<word>("name", "EnSight");
    if (!ensightDir.isAbsolute())
    {
        ensightDir = args.rootPath()/args.globalCaseName()/ensightDir;
    }


    //
    // Output configuration (geometry related)
    //
    ensightMesh::options writeOpts(format);
    writeOpts.noPatches(args.optionFound("noPatches"));

    if (args.optionFound("patches"))
    {
        writeOpts.patchSelection(args.optionReadList<wordRe>("patches"));
    }
    if (args.optionFound("faceZones"))
    {
        writeOpts.faceZoneSelection(args.optionReadList<wordRe>("faceZones"));
    }

    //
    // output configuration (field related)
    //
    const bool noLagrangian = args.optionFound("noLagrangian");

    wordReList fieldPatterns;
    if (args.optionFound("fields"))
    {
        fieldPatterns = args.optionReadList<wordRe>("fields");
    }

    word cellZoneName;
    if (args.optionReadIfPresent("cellZone", cellZoneName))
    {
        Info<< "Converting cellZone " << cellZoneName
            << " only (puts outside faces into patch "
            << mesh.boundaryMesh()[0].name() << ")"
            << endl;
    }
    meshSubsetHelper meshRef(mesh, meshSubsetHelper::ZONE, cellZoneName);

    //
    // Open new ensight case file, initialize header etc.
    //
    ensightCase ensCase
    (
        ensightDir,
        args.globalCaseName(),
        caseOpts
    );


    // Construct the Ensight mesh
    ensightMesh ensMesh(meshRef.mesh(), writeOpts);

    if (Pstream::master())
    {
        Info<< "Converting " << timeDirs.size() << " time steps" << nl;
        ensCase.printInfo(Info) << endl;
    }

    #include "checkMeshMoving.H"
    #include "findCloudFields.H"

    // test the pre-check variable if there is a moving mesh
    // time-set for geometries
    // TODO: split off into separate time-set,
    // but need to verify ensight spec

    Info<< "Startup in "
        << timer.cpuTimeIncrement() << " s, "
        << mem.update().size() << " kB" << nl << endl;

    // Get the list of supported classes/fields
    HashTable<wordHashSet> usableObjects;
    {
        // Initially all possible objects that are available at the final time
        IOobjectList objects(mesh, timeDirs.last().name());

        // Categorize by classes, pre-filter on name (if requested)
        usableObjects =
        (
            fieldPatterns.empty()
          ? objects.classes()
          : objects.classes(fieldPatterns)
        );

        // Limit to types that we explicitly handle
        usableObjects.filterKeys(volFieldTypes);

        // Force each field-type into existence (simplifies code logic
        // and doesn't cost much) and simultaneously remove all
        // "*_0" restart fields

        for (auto fieldType : volFieldTypes)
        {
            usableObjects
            (
                fieldType
            ).filterKeys
            (
                [](const word& k){ return k.endsWith("_0"); },
                true  // prune
            );
        }
    }

    // ignore special fields (_0 fields),
    // ignore fields we don't handle,
    // ignore fields that are not available for all time-steps
    HashTable<bool> fieldsToUse;

    forAll(timeDirs, timeIndex)
    {
        runTime.setTime(timeDirs[timeIndex], timeIndex);
        ensCase.nextTime(timeDirs[timeIndex]);

        Info<< "Time [" << timeIndex << "] = " << runTime.timeName() << nl;

        polyMesh::readUpdateState meshState = mesh.readUpdate();
        if (meshState != polyMesh::UNCHANGED)
        {
            meshRef.correct();
            ensMesh.expire();
            ensMesh.correct();
        }

        if (timeIndex == 0 || meshMoving)
        {
            autoPtr<ensightGeoFile> os = ensCase.newGeometry(meshMoving);
            ensMesh.write(os);
        }


        // Cell field data output
        // ~~~~~~~~~~~~~~~~~~~~~~
        Info<< "Write volume field (";

        for (auto fieldType : volFieldTypes)
        {
            // For convenience, just force each field-type into existence.
            // This simplifies code logic and doesn't cost much at all.
            wordHashSet& fieldNames = usableObjects(fieldType);

            forAllIters(fieldNames, fieldIter)
            {
                const word& fieldName = fieldIter.key();

                #include "checkData.H"

                // Partially complete field?
                if (!fieldsToUse[fieldName])
                {
                    fieldNames.erase(fieldIter);
                    continue;
                }

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
                    wrote = ensightOutput::writeField<scalar>
                    (
                        meshRef.interpolate(vf),
                        ensMesh,
                        os,
                        nodeValues
                    );
                }
                else if (fieldType == volVectorField::typeName)
                {
                    autoPtr<ensightFile> os = ensCase.newData<vector>
                    (
                        fieldName
                    );

                    volVectorField vf(fieldObject, mesh);
                    wrote = ensightOutput::writeField<vector>
                    (
                        meshRef.interpolate(vf),
                        ensMesh,
                        os,
                        nodeValues
                    );
                }
                else if (fieldType == volSphericalTensorField::typeName)
                {
                    autoPtr<ensightFile> os = ensCase.newData<sphericalTensor>
                    (
                        fieldObject.name()
                    );

                    volSphericalTensorField vf(fieldObject, mesh);
                    wrote = ensightOutput::writeField<sphericalTensor>
                    (
                        meshRef.interpolate(vf),
                        ensMesh,
                        os,
                        nodeValues
                    );
                }
                else if (fieldType == volSymmTensorField::typeName)
                {
                    autoPtr<ensightFile> os = ensCase.newData<symmTensor>
                    (
                        fieldName
                    );

                    volSymmTensorField vf(fieldObject, mesh);
                    wrote = ensightOutput::writeField<symmTensor>
                    (
                        meshRef.interpolate(vf),
                        ensMesh,
                        os,
                        nodeValues
                    );
                }
                else if (fieldType == volTensorField::typeName)
                {
                    autoPtr<ensightFile> os = ensCase.newData<tensor>
                    (
                        fieldName
                    );

                    volTensorField vf(fieldObject, mesh);
                    wrote = ensightOutput::writeField<tensor>
                    (
                        meshRef.interpolate(vf),
                        ensMesh,
                        os,
                        nodeValues
                    );
                }
                // DimensionedFields
                else if
                (
                    fieldType
                 == volScalarField::Internal::typeName
                )
                {
                    autoPtr<ensightFile> os = ensCase.newData<scalar>
                    (
                        fieldName
                    );

                    volScalarField::Internal df
                    (
                        fieldObject,
                        mesh
                    );
                    wrote = ensightOutput::writeField<scalar>
                    (
                        meshRef.interpolate<scalar>(df),
                        ensMesh,
                        os,
                        nodeValues
                    );
                }
                else if
                (
                    fieldType
                 == volVectorField::Internal::typeName
                )
                {
                    autoPtr<ensightFile> os = ensCase.newData<vector>
                    (
                        fieldName
                    );

                    volVectorField::Internal df
                    (
                        fieldObject,
                        mesh
                    );
                    wrote = ensightOutput::writeField<vector>
                    (
                        meshRef.interpolate<vector>(df),
                        ensMesh,
                        os,
                        nodeValues
                    );
                }
                else if
                (
                    fieldType
                 == volSphericalTensorField::Internal::typeName
                )
                {
                    autoPtr<ensightFile> os = ensCase.newData<sphericalTensor>
                    (
                        fieldName
                    );

                    volSphericalTensorField::Internal df
                    (
                        fieldObject,
                        mesh
                    );
                    wrote = ensightOutput::writeField<sphericalTensor>
                    (
                        meshRef.interpolate<sphericalTensor>(df),
                        ensMesh,
                        os,
                        nodeValues
                    );
                }
                else if
                (
                    fieldType
                 == volSymmTensorField::Internal::typeName
                )
                {
                    autoPtr<ensightFile> os = ensCase.newData<symmTensor>
                    (
                        fieldName
                    );

                    volSymmTensorField::Internal df
                    (
                        fieldObject,
                        mesh
                    );
                    wrote = ensightOutput::writeField<symmTensor>
                    (
                        meshRef.interpolate<symmTensor>(df),
                        ensMesh,
                        os,
                        nodeValues
                    );
                }
                else if
                (
                    fieldType
                 == volTensorField::Internal::typeName
                )
                {
                    autoPtr<ensightFile> os = ensCase.newData<tensor>
                    (
                        fieldName
                    );

                    volTensorField::Internal df
                    (
                        fieldObject,
                        mesh
                    );
                    wrote = ensightOutput::writeField<tensor>
                    (
                        meshRef.interpolate<tensor>(df),
                        ensMesh,
                        os,
                        nodeValues
                    );
                }
                else
                {
                    // Do not currently handle this type
                    // - blacklist for the future.
                    fieldsToUse.set(fieldName, false);
                }

                if (wrote)
                {
                    Info<< ' ' << fieldName;
                }
            }
        }
        Info<< " )" << nl;


        // Cloud field data output
        // ~~~~~~~~~~~~~~~~~~~~~~~

        forAll(cloudNames, cloudNo)
        {
            const word& cloudName = cloudNames[cloudNo];
            const HashTable<word>& theseCloudFields = cloudFields[cloudName];

            fileNameList currentCloudDirs = readDir
            (
                runTime.timePath()/regionPrefix/cloud::prefix,
                fileName::DIRECTORY
            );

            Info<< "Write " << cloudName << " (";

            bool cloudExists = currentCloudDirs.found(cloudName);
            reduce(cloudExists, orOp<bool>());

            {
                autoPtr<ensightFile> os = ensCase.newCloud(cloudName);

                ensightCloud::writePositions
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

            forAllConstIter(HashTable<word>, theseCloudFields, fieldIter)
            {
                const word& fieldName = fieldIter.key();
                const word& fieldType = fieldIter();

                IOobject fieldObject
                (
                    fieldName,
                    mesh.time().timeName(),
                    cloud::prefix/cloudName,
                    mesh,
                    IOobject::MUST_READ
                );

                // cannot have field without cloud positions
                bool fieldExists = cloudExists;
                if (cloudExists)
                {
                    fieldExists =
                        fieldObject.typeHeaderOk<IOField<scalar>>(false);

                    reduce(fieldExists, orOp<bool>());
                }

                bool wrote = false;
                if (fieldType == scalarIOField::typeName)
                {
                    autoPtr<ensightFile> os =
                        ensCase.newCloudData<scalar>(cloudName, fieldName);

                    wrote = ensightCloud::writeCloudField<scalar>
                    (
                        fieldObject, fieldExists, os
                    );
                }
                else if (fieldType == vectorIOField::typeName)
                {
                    autoPtr<ensightFile> os =
                        ensCase.newCloudData<vector>(cloudName, fieldName);

                    wrote = ensightCloud::writeCloudField<vector>
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

        Info<< "Wrote in "
            << timer.cpuTimeIncrement() << " s, "
            << mem.update().size() << " kB" << nl << nl;
    }

    ensCase.write();

    Info<< "End: "
        << timer.elapsedCpuTime() << " s, "
        << mem.update().peak() << " kB (peak)" << nl << endl;

    return 0;
}


// ************************************************************************* //
