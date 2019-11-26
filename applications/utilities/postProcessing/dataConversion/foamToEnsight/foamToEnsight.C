/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2018 OpenCFD Ltd.
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
    Translate OpenFOAM data to EnSight format.
    An Ensight part is created for the internalMesh and for each patch.

Usage
    \b foamToEnsight [OPTION]

    Options:
      - \par -ascii
        Write Ensight data in ASCII format instead of "C Binary"

      - \par -fields \<fields\>
        Specify single or multiple fields to write (all by default)
        For example,
        \verbatim
          -fields T
          -fields '(p T U \"alpha.*\")'
        \endverbatim
        The quoting is required to avoid shell expansions and to pass the
        information as a single argument.

      - \par -no-boundary
        Suppress writing any patches.

      - \par -no-internal
        Suppress writing the internal mesh.

      - \par -no-lagrangian
        Suppress writing lagrangian positions and fields.

      - \par -patches patch or patch list
        Specify particular patches to write.

      - \par -faceZones zone or zone list
        Specify faceZones to write, with wildcards

      - \par -cellZone zoneName
        Specify single cellZone to write (not lagrangian)

      - \par -noZero
        Exclude the often incomplete initial conditions.

      - \par -name \<subdir\>
        Define sub-directory name to use for Ensight data (default: "EnSight")

      - \par -width \<n\>
        Width of Ensight data subdir (default: 8)

Note
    Writes to \a EnSight directory to avoid collisions with
    foamToEnsightParts

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "IOobjectList.H"
#include "IOmanip.H"
#include "OFstream.H"
#include "PstreamCombineReduceOps.H"
#include "HashOps.H"

#include "fvc.H"
#include "fieldTypes.H"
#include "volFields.H"
#include "scalarIOField.H"
#include "vectorIOField.H"

// file-format/conversion
#include "ensightCase.H"
#include "ensightGeoFile.H"
#include "ensightMesh.H"
#include "ensightOutputCloud.H"
#include "ensightOutputVolField.H"
#include "fvMeshSubsetProxy.H"

// local files
#include "readFields.H"
#include "writeVolFields.H"
#include "writeDimFields.H"

#include "memInfo.H"

using namespace Foam;

//- Get internal field and make it a zero-gradient volume field with subsetting
template<class GeoField>
tmp<GeoField>
getZeroGradInternalField(IOobject& io, const fvMeshSubsetProxy& proxy)
{
    auto tfield = tmp<typename GeoField::Internal>::New(io, proxy.baseMesh());
    return proxy.interpolateInternal(tfield);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Translate OpenFOAM data to Ensight format with a part for"
        " the internalMesh and for each patch."
    );
    timeSelector::addOptions();

    // Less frequently used - reduce some clutter
    argList::setAdvanced("decomposeParDict");
    argList::setAdvanced("noFunctionObjects");

    #include "addRegionOption.H"

    argList::addBoolOption
    (
        "ascii",
        "Write in ASCII format instead of 'C Binary'"
    );
    argList::addBoolOption
    (
        "nodeValues",
        "Write values at nodes"
    );
    argList::addBoolOption
    (
        "no-boundary",  // noPatches
        "Suppress writing any patches"
    );
    argList::addOptionCompat("no-boundary", {"noPatches", 1806});

    argList::addBoolOption
    (
        "no-internal",
        "Suppress writing the internal mesh"
    );
    argList::addBoolOption
    (
        "no-lagrangian",  // noLagrangian
        "Suppress writing lagrangian positions and fields"
    );
    argList::addOptionCompat("no-lagrangian", {"noLagrangian", 1806});

    argList::addOption
    (
        "patches",
        "wordRes",
        "Specify single patch or multiple patches to write\n"
        "Eg, 'inlet' or '(outlet \"inlet.*\")'"
    );
    argList::addOption
    (
        "faceZones",
        "wordRes",
        "Specify single or multiple faceZones to write\n"
        "Eg, 'cells' or '( slice \"mfp-.*\" )'."
    );
    argList::addOption
    (
        "fields",
        "wordRes",
        "Specify single or multiple fields to write (all by default)\n"
        "Eg, 'T' or '( \"U.*\" )'"
    );
    argList::addOption
    (
        "cellZone",
        "word",
        "Specify cellZone to write"
    );
    argList::addOption
    (
        "name",
        "subdir",
        "Sub-directory name for Ensight output (default: 'EnSight')"
    );
    argList::addOption
    (
        "width",
        "n",
        "Width of Ensight data subdir"
    );

    #include "setRootCase.H"

    // Default to binary output, unless otherwise specified
    const IOstream::streamFormat format =
    (
        args.found("ascii")
      ? IOstream::ASCII
      : IOstream::BINARY
    );

    const bool nodeValues = args.found("nodeValues");

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

    caseOpts.nodeValues(args.found("nodeValues"));
    caseOpts.width(args.get<label>("width", 8));
    caseOpts.overwrite(true); // remove existing output directory

    // Can also have separate directory for lagrangian
    // caseOpts.separateCloud(true);


    // Define sub-directory name to use for EnSight data.
    // The path to the ensight directory is at case level only
    // - For parallel cases, data only written from master
    fileName outputDir = args.get<word>("name", "EnSight");
    if (!outputDir.isAbsolute())
    {
        outputDir = args.globalPath()/outputDir;
    }


    //
    // Output configuration (geometry related)
    //
    ensightMesh::options writeOpts(format);
    writeOpts.useBoundaryMesh(!args.found("no-boundary"));
    writeOpts.useInternalMesh(!args.found("no-internal"));
    const bool doLagrangian = !args.found("no-lagrangian");

    if (args.found("patches"))
    {
        writeOpts.patchSelection(args.getList<wordRe>("patches"));
    }
    if (args.found("faceZones"))
    {
        writeOpts.faceZoneSelection(args.getList<wordRe>("faceZones"));
    }

    //
    // Output configuration (field related)
    //

    wordRes fieldPatterns;
    args.readListIfPresent<wordRe>("fields", fieldPatterns);

    word cellZoneName;
    if (args.readIfPresent("cellZone", cellZoneName))
    {
        Info<< "Converting cellZone " << cellZoneName
            << " only, with new outside faces as \"oldInternalFaces\"."
            << nl;
    }

    // Ignored (unproxied) if cellZoneName is empty
    fvMeshSubsetProxy meshProxy(mesh, fvMeshSubsetProxy::ZONE, cellZoneName);

    // New ensight case file, initialize header etc.
    ensightCase ensCase(outputDir, args.globalCaseName(), caseOpts);

    // Construct the Ensight mesh
    ensightMesh ensMesh(meshProxy.mesh(), writeOpts);

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


    // Initially all possible objects that are available at the final time
    wordHashSet testedObjectNames;
    {
        IOobjectList objects(mesh, timeDirs.last().name());

        if (!fieldPatterns.empty())
        {
            objects.filterObjects(fieldPatterns);
        }

        // Remove "*_0" restart fields
        objects.prune_0();

        // Only retain volume and dimensioned fields.
        objects.filterClasses
        (
            [](const word& clsName){
                return
                (
                    fieldTypes::volume.found(clsName)
                 || fieldTypes::internal.found(clsName)
                );
            }
        );

        wordList objectNames(objects.sortedNames());

        // Check availability for all times...
        checkData(meshProxy.baseMesh(), timeDirs, objectNames);

        testedObjectNames = objectNames;
    }


    forAll(timeDirs, timeIndex)
    {
        runTime.setTime(timeDirs[timeIndex], timeIndex);
        ensCase.nextTime(timeDirs[timeIndex]);

        Info<< "Time [" << timeIndex << "] = " << runTime.timeName() << nl;

        polyMesh::readUpdateState meshState = mesh.readUpdate();
        if (meshState != polyMesh::UNCHANGED)
        {
            meshProxy.correct();
            ensMesh.expire();
            ensMesh.correct();
        }

        if (timeIndex == 0 || meshMoving)
        {
            autoPtr<ensightGeoFile> os = ensCase.newGeometry(meshMoving);
            ensMesh.write(os);
        }

        // Objects at this time
        IOobjectList objects(meshProxy.baseMesh(), runTime.timeName());

        // Restrict to objects that are available for all times
        objects.filterObjects(testedObjectNames);

        // Volume, internal, point fields
        #include "convertVolumeFields.H"

        // Write lagrangian data
        #include "convertLagrangian.H"

        Info<< "Wrote in "
            << timer.cpuTimeIncrement() << " s, "
            << mem.update().size() << " kB" << nl << nl;
    }

    ensCase.write();

    Info<< "\nEnd: "
        << timer.elapsedCpuTime() << " s, "
        << mem.update().peak() << " kB (peak)" << nl << endl;

    return 0;
}


// ************************************************************************* //
