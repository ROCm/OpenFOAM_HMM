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
    foamToEnsight

Group
    grpPostProcessingUtilitie

Description
    Translates OpenFOAM data to EnSight format.

    An Ensight part is created for the internalMesh and for each patch.

Usage
    - foamToEnsight [OPTION] \n
    Translates OpenFOAM data to EnSight format

    \param -ascii \n
    Write Ensight data in ASCII format instead of "C Binary"

    \param -noZero \n
    Exclude the often incomplete initial conditions.

    \param -patches patchList \n
    Specify particular patches to write.
    Specifying an empty list suppresses writing the internalMesh.

    \param -noPatches \n
    Suppress writing any patches.

    \param -faceZones zoneList \n
    Specify faceZones to write, with wildcards

    \param -cellZone zoneName \n
    Specify single cellZone to write (not lagrangian)

    \param -width \<n\>\n
    Width of EnSight data subdir (default: 8)

Note
    Parallel support for cloud data is not supported
    - writes to \a EnSight directory to avoid collisions with foamToEnsightParts

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "IOobjectList.H"
#include "IOmanip.H"
#include "OFstream.H"

#include "volFields.H"

#include "labelIOField.H"
#include "scalarIOField.H"
#include "tensorIOField.H"

#include "ensightFile.H"
#include "ensightMesh.H"
#include "ensightField.H"
#include "ensightCloud.H"

#include "fvc.H"
#include "cellSet.H"
#include "fvMeshSubset.H"

#include "memInfo.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool inFileNameList
(
    const fileNameList& nameList,
    const word& name
)
{
    forAll(nameList, i)
    {
        if (nameList[i] == name)
        {
            return true;
        }
    }

    return false;
}



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
        "define sub-directory name to use for ensight data "
        "(default: 'EnSight')"
    );
    argList::addOption
    (
        "width",
        "n",
        "width of ensight data subdir"
    );

    // the volume field types that we handle
    const label nVolFieldTypes = 10;
    const word volFieldTypes[] =
    {
        volScalarField::typeName,
        volVectorField::typeName,
        volSphericalTensorField::typeName,
        volSymmTensorField::typeName,
        volTensorField::typeName,

        volScalarField::DimensionedInternalField::typeName,
        volVectorField::DimensionedInternalField::typeName,
        volSphericalTensorField::DimensionedInternalField::typeName,
        volSymmTensorField::DimensionedInternalField::typeName,
        volTensorField::DimensionedInternalField::typeName
    };

    #include "setRootCase.H"

    // default to binary output, unless otherwise specified
    const IOstream::streamFormat format =
    (
        args.optionFound("ascii")
      ? IOstream::ASCII
      : IOstream::BINARY
    );

    const bool nodeValues = args.optionFound("nodeValues");

    cpuTime timer;
    memInfo mem;
    Info<< "Initial memory "
        << mem.update().size() << " kB" << endl;

    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    // adjust output width
    if (args.optionFound("width"))
    {
        ensightFile::subDirWidth(args.optionRead<label>("width"));
    }

    // define sub-directory name to use for EnSight data
    fileName ensightDir = "EnSight";
    args.optionReadIfPresent("name", ensightDir);

    // Path to EnSight directory at case level only
    // - For parallel cases, data only written from master
    if (!ensightDir.isAbsolute())
    {
        ensightDir = args.rootPath()/args.globalCaseName()/ensightDir;
    }

    const fileName dataDir  = ensightDir/"data";
    const fileName dataMask = dataDir.name()/ensightFile::mask();

    if (Pstream::master())
    {
        // EnSight and EnSight/data directories must exist
        // - remove old data for a clean conversion of everything
        if (isDir(ensightDir))
        {
            rmDir(ensightDir);
        }

        mkDir(dataDir);
    }

    #include "createNamedMesh.H"

    // Mesh instance (region0 gets filtered out)
    fileName regionPrefix;
    if (regionName != polyMesh::defaultRegion)
    {
        regionPrefix = regionName;
    }

    // Start of case file header output
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    OFstream *ensightCaseFilePtr(nullptr);
    if (Pstream::master())
    {
        fileName caseFileName = args.globalCaseName() + ".case";

        Info<< "Converting " << timeDirs.size() << " time steps" << nl
            << "Ensight case: " << caseFileName.c_str() << endl;

        // The case file is always ASCII
        ensightCaseFilePtr = new OFstream
        (
            ensightDir/caseFileName,
            IOstream::ASCII
        );

        ensightCaseFilePtr->setf(ios_base::left);
        ensightCaseFilePtr->setf(ios_base::scientific, ios_base::floatfield);
        ensightCaseFilePtr->precision(5);

        *ensightCaseFilePtr
            << "FORMAT" << nl
            << "type: ensight gold" << nl << nl;
    }

    OFstream& ensightCaseFile = *ensightCaseFilePtr;

    // Construct the EnSight mesh
    const bool selectedPatches = args.optionFound("patches");
    wordReList patchPatterns;
    if (selectedPatches)
    {
        patchPatterns = wordReList(args.optionLookup("patches")());
    }
    const bool selectedZones = args.optionFound("faceZones");
    wordReList zonePatterns;
    if (selectedZones)
    {
        zonePatterns = wordReList(args.optionLookup("faceZones")());
    }

    const bool selectedFields = args.optionFound("fields");
    wordReList fieldPatterns;
    if (selectedFields)
    {
        fieldPatterns = wordReList(args.optionLookup("fields")());
    }

    word cellZoneName;
    const bool doCellZone = args.optionReadIfPresent("cellZone", cellZoneName);

    fvMeshSubset meshSubsetter(mesh);
    if (doCellZone)
    {
        Info<< "Converting cellZone " << cellZoneName
            << " only (puts outside faces into patch "
            << mesh.boundaryMesh()[0].name()
            << ")" << endl;
        const cellZone& cz = mesh.cellZones()[cellZoneName];
        cellSet c0(mesh, "c0", labelHashSet(cz));
        meshSubsetter.setLargeCellSubset(c0, 0);
    }

    ensightMesh eMesh
    (
        (
            meshSubsetter.hasSubMesh()
          ? meshSubsetter.subMesh()
          : meshSubsetter.baseMesh()
        ),
        args.optionFound("noPatches"),
        selectedPatches,
        patchPatterns,
        selectedZones,
        zonePatterns,
        format
    );

    // Set Time to the last time before looking for the lagrangian objects
    runTime.setTime(timeDirs.last(), timeDirs.size()-1);

    IOobjectList objects(mesh, runTime.timeName());

    #include "checkMeshMoving.H"

    if (Pstream::master())
    {
        // test the pre-check variable if there is a moving mesh
        // time-set for geometries
        // TODO: split off into separate time-set,
        // but need to verify ensight spec

        if (meshMoving)
        {
            ensightCaseFile
                << "GEOMETRY" << nl
                << setw(16) << "model: 1"
                << (dataMask/ensightMesh::geometryName).c_str() << nl;
        }
        else
        {
            ensightCaseFile
                << "GEOMETRY" << nl
                << setw(16) << "model:"
                << ensightMesh::geometryName << nl;
        }
    }

    // Identify if lagrangian data exists any time step, and add clouds
    // to the 'cloudNames' (sorted list)
    wordList cloudNames;
    {
        wordHashSet allCloudNames;

        forAll(timeDirs, timeI)
        {
            runTime.setTime(timeDirs[timeI], timeI);

            fileNameList cloudDirs = readDir
            (
                runTime.timePath()/regionPrefix/cloud::prefix,
                fileName::DIRECTORY
            );

            forAll(cloudDirs, cloudI)
            {
                IOobjectList cloudObjs
                (
                    mesh,
                    runTime.timeName(),
                    cloud::prefix/cloudDirs[cloudI]
                );

                IOobject* positionsPtr = cloudObjs.lookup(word("positions"));

                if (positionsPtr)
                {
                    allCloudNames.insert(cloudDirs[cloudI]);
                }
            }
        }

        // sorted for consistency
        cloudNames = allCloudNames.sortedToc();
    }

    // ignore special fields (_0 fields),
    // ignore fields we don't handle,
    // ignore fields that are not available for all time-steps
    HashTable<bool> fieldsToUse;

    HashTable<HashTable<word>> allCloudFields;
    forAll(cloudNames, cloudNo)
    {
        const word& cloudName = cloudNames[cloudNo];

        // Add the name of the cloud(s) to the case file header
        if (Pstream::master())
        {
            ensightCaseFile
                << setw(16) << "measured: 1"
                << fileName(dataMask/cloud::prefix/cloudName/"positions").c_str()
                << nl;
        }

        // Create a new hash table for each cloud
        allCloudFields.insert(cloudName, HashTable<word>());

        // Identify the new cloud in the hash table
        HashTable<HashTable<word>>::iterator newCloudIter =
            allCloudFields.find(cloudName);

        // Loop over all times to build list of fields and field types
        // for each cloud
        forAll(timeDirs, timeI)
        {
            runTime.setTime(timeDirs[timeI], timeI);

            IOobjectList cloudObjs
            (
                mesh,
                runTime.timeName(),
                cloud::prefix/cloudName
            );

            forAllConstIter(IOobjectList, cloudObjs, fieldIter)
            {
                const IOobject obj = *fieldIter();

                if (obj.name() != "positions")
                {
                    // Add field and field type
                    newCloudIter().insert
                    (
                        obj.name(),
                        obj.headerClassName()
                    );
                }
            }
        }
    }

    Info<< "Startup in "
        << timer.cpuTimeIncrement() << " s, "
        << mem.update().size() << " kB" << nl << endl;

    label nTimeSteps = 0;
    forAll(timeDirs, timeIndex)
    {
        ++nTimeSteps;
        runTime.setTime(timeDirs[timeIndex], timeIndex);

        Info<< "Time [" << timeIndex << "] = " << runTime.timeName() << nl;

        if (Pstream::master())
        {
            // the data/ITER subdirectory must exist
            // Note that data/ITER is indeed a valid ensight::FileName
            const fileName subDir = ensightFile::subDir(timeIndex);
            mkDir(dataDir/subDir);

            // place a timestamp in the directory for future reference
            OFstream timeStamp(dataDir/subDir/"time");
            timeStamp
                << "#   timestep time" << nl
                << subDir.c_str() << " " << runTime.timeName() << nl;
        }

        polyMesh::readUpdateState meshState = mesh.readUpdate();
        if (timeIndex != 0 && meshSubsetter.hasSubMesh())
        {
            Info<< "Converting cellZone " << cellZoneName
                << " only (puts outside faces into patch "
                << mesh.boundaryMesh()[0].name()
                << ")" << endl;
            const cellZone& cz = mesh.cellZones()[cellZoneName];
            cellSet c0(mesh, "c0", labelHashSet(cz));
            meshSubsetter.setLargeCellSubset(c0, 0);
        }

        if (meshState != polyMesh::UNCHANGED)
        {
            eMesh.correct();
        }

        if (timeIndex == 0 || meshMoving)
        {
            eMesh.write
            (
                dataDir,
                timeIndex,
                meshMoving,
                ensightCaseFile
            );
        }


        // Start of field data output
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (timeIndex == 0 && Pstream::master())
        {
            ensightCaseFile<< nl << "VARIABLE" << nl;
        }


        // Cell field data output
        // ~~~~~~~~~~~~~~~~~~~~~~
        Info<< "Write volume field (";

        for (label i=0; i<nVolFieldTypes; ++i)
        {
            wordList fieldNames = objects.names(volFieldTypes[i]);

            forAll(fieldNames, j)
            {
                const word& fieldName = fieldNames[j];

                // Check if the field has to be exported
                if (selectedFields)
                {
                    if (!findStrings(fieldPatterns, fieldName))
                    {
                        continue;
                    }
                }

                #include "checkData.H"

                if (!fieldsToUse[fieldName])
                {
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

                if (volFieldTypes[i] == volScalarField::typeName)
                {
                    volScalarField vf(fieldObject, mesh);
                    ensightField<scalar>
                    (
                        volField(meshSubsetter, vf),
                        eMesh,
                        dataDir,
                        timeIndex,
                        nodeValues,
                        ensightCaseFile
                    );
                }
                else if (volFieldTypes[i] == volVectorField::typeName)
                {
                    volVectorField vf(fieldObject, mesh);
                    ensightField<vector>
                    (
                        volField(meshSubsetter, vf),
                        eMesh,
                        dataDir,
                        timeIndex,
                        nodeValues,
                        ensightCaseFile
                    );
                }
                else if (volFieldTypes[i] == volSphericalTensorField::typeName)
                {
                    volSphericalTensorField vf(fieldObject, mesh);
                    ensightField<sphericalTensor>
                    (
                        volField(meshSubsetter, vf),
                        eMesh,
                        dataDir,
                        timeIndex,
                        nodeValues,
                        ensightCaseFile
                    );
                }
                else if (volFieldTypes[i] == volSymmTensorField::typeName)
                {
                    volSymmTensorField vf(fieldObject, mesh);
                    ensightField<symmTensor>
                    (
                        volField(meshSubsetter, vf),
                        eMesh,
                        dataDir,
                        timeIndex,
                        nodeValues,
                        ensightCaseFile
                    );
                }
                else if (volFieldTypes[i] == volTensorField::typeName)
                {
                    volTensorField vf(fieldObject, mesh);
                    ensightField<tensor>
                    (
                        volField(meshSubsetter, vf),
                        eMesh,
                        dataDir,
                        timeIndex,
                        nodeValues,
                        ensightCaseFile
                    );
                }
                // DimensionedFields
                else if
                (
                    volFieldTypes[i]
                 == volScalarField::DimensionedInternalField::typeName
                )
                {
                    volScalarField::DimensionedInternalField df
                    (
                        fieldObject,
                        mesh
                    );
                    ensightField<scalar>
                    (
                        volField<scalar>(meshSubsetter, df),
                        eMesh,
                        dataDir,
                        timeIndex,
                        nodeValues,
                        ensightCaseFile
                    );
                }
                else if
                (
                    volFieldTypes[i]
                 == volVectorField::DimensionedInternalField::typeName
                )
                {
                    volVectorField::DimensionedInternalField df
                    (
                        fieldObject,
                        mesh
                    );
                    ensightField<vector>
                    (
                        volField<vector>(meshSubsetter, df),
                        eMesh,
                        dataDir,
                        timeIndex,
                        nodeValues,
                        ensightCaseFile
                    );
                }
                else if
                (
                    volFieldTypes[i]
                 == volSphericalTensorField::DimensionedInternalField::typeName
                )
                {
                    volSphericalTensorField::DimensionedInternalField df
                    (
                        fieldObject,
                        mesh
                    );
                    ensightField<sphericalTensor>
                    (
                        volField<sphericalTensor>(meshSubsetter, df),
                        eMesh,
                        dataDir,
                        timeIndex,
                        nodeValues,
                        ensightCaseFile
                    );
                }
                else if
                (
                    volFieldTypes[i]
                 == volSymmTensorField::DimensionedInternalField::typeName
                )
                {
                    volSymmTensorField::DimensionedInternalField df
                    (
                        fieldObject,
                        mesh
                    );
                    ensightField<symmTensor>
                    (
                        volField<symmTensor>(meshSubsetter, df),
                        eMesh,
                        dataDir,
                        timeIndex,
                        nodeValues,
                        ensightCaseFile
                    );
                }
                else if
                (
                    volFieldTypes[i]
                 == volTensorField::DimensionedInternalField::typeName
                )
                {
                    volTensorField::DimensionedInternalField df
                    (
                        fieldObject,
                        mesh
                    );
                    ensightField<tensor>
                    (
                        volField<tensor>(meshSubsetter, df),
                        eMesh,
                        dataDir,
                        timeIndex,
                        nodeValues,
                        ensightCaseFile
                    );
                }
                else
                {
                    // Do not currently handle this type - blacklist for the future.
                    fieldsToUse.set(fieldName, false);
                }
            }
        }
        Info<< " )" << nl;


        // Cloud field data output
        // ~~~~~~~~~~~~~~~~~~~~~~~

        forAll(cloudNames, cloudNo)
        {
            const word& cloudName = cloudNames[cloudNo];
            if (!allCloudFields.found(cloudName))
            {
                // extra safety
                continue;
            }

            const HashTable<word>& cloudFields = allCloudFields[cloudName];

            fileNameList currentCloudDirs = readDir
            (
                runTime.timePath()/regionPrefix/cloud::prefix,
                fileName::DIRECTORY
            );

            Info<< "Write " << cloudName << " (";

            bool cloudExists = inFileNameList(currentCloudDirs, cloudName);
            reduce(cloudExists, orOp<bool>());

            ensightParticlePositions
            (
                mesh,
                dataDir,
                timeIndex,
                cloudName,
                cloudExists,
                format
            );

            forAllConstIter(HashTable<word>, cloudFields, fieldIter)
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

                bool fieldExists = fieldObject.typeHeaderOk<IOField<scalar>>
                (
                    false
                );
                reduce(fieldExists, orOp<bool>());

                if (fieldType == scalarIOField::typeName)
                {
                    ensightCloudField<scalar>
                    (
                        fieldObject,
                        dataDir,
                        timeIndex,
                        cloudName,
                        cloudNo,
                        ensightCaseFile,
                        fieldExists,
                        format
                    );
                }
                else if (fieldType == vectorIOField::typeName)
                {
                    ensightCloudField<vector>
                    (
                        fieldObject,
                        dataDir,
                        timeIndex,
                        cloudName,
                        cloudNo,
                        ensightCaseFile,
                        fieldExists,
                        format
                    );
                }
            }
            Info<< " )" << nl;
        }

        Info<< "Wrote in "
            << timer.cpuTimeIncrement() << " s, "
            << mem.update().size() << " kB" << nl << nl;
    }

    #include "ensightCaseTail.H"

    if (ensightCaseFilePtr) // on master only
    {
        delete ensightCaseFilePtr;
    }

    Info<< "End: "
        << timer.elapsedCpuTime() << " s, "
        << mem.update().peak() << " kB (peak)" << nl << endl;

    return 0;
}


// ************************************************************************* //
