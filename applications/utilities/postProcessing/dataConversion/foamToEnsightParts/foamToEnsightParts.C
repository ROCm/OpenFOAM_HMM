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
    foamToEnsightParts

Group
    grpPostProcessingUtilities

Description
    Translate OpenFOAM data to Ensight format with an Ensight part
    for each cellZone and patch.

Usage
    \b foamToEnsightParts [OPTION]

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

      - \par -noZero
        Exclude the often incomplete initial conditions.

      - \par -index \<start\>
        Ignore the time index contained in the time file and use a
        simple indexing when creating the \c Ensight/data/######## files.

      - \par -no-lagrangian
        Suppress writing lagrangian positions and fields.

      - \par -no-mesh
        Suppress writing the geometry. Can be useful for converting partial
        results for a static geometry.

      - \par -noZero
        Exclude the often incomplete initial conditions.

      - \par -name \<subdir\>
        Define sub-directory name to use for Ensight data (default: "Ensight")

      - \par -width \<n\>
        Width of Ensight data subdir

Note
    - no parallel data.
    - writes to \a Ensight directory to avoid collisions with foamToEnsight.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "IOobjectList.H"
#include "IOmanip.H"
#include "OFstream.H"
#include "PstreamCombineReduceOps.H"
#include "HashOps.H"

#include "fieldTypes.H"
#include "volFields.H"
#include "scalarIOField.H"
#include "vectorIOField.H"

// file-format/conversion
#include "ensightCase.H"
#include "ensightGeoFile.H"
#include "ensightParts.H"
#include "ensightOutputCloud.H"
#include "ensightOutputVolField.H"
#include "fvMeshSubsetProxy.H"

// local files
#include "readFields.H"
#include "writeVolFields.H"
#include "writeDimFields.H"

#include "memInfo.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Translate OpenFOAM data to Ensight format with an Ensight part"
        " for each cellZone and patch."
    );

    // Enable -constant
    // Probably don't need -withZero though, since the fields are vetted
    // afterwards anyhow
    timeSelector::addOptions(true, false); // constant(true), zero(false)
    #include "addRegionOption.H"

    argList::noParallel();
    argList::addBoolOption
    (
        "ascii",
        "Write in ASCII format instead of 'C Binary'"
    );
    argList::addOption
    (
        "index",
        "start",
        "Ignore the time index contained in the uniform/time file"
        " and use simple indexing when creating files"
    );
    argList::addBoolOption
    (
        "no-lagrangian", // noLagrangian
        "Suppress writing lagrangian positions and fields"
    );
    argList::addOptionCompat("no-lagrangian", {"noLagrangian", 1806});

    argList::addBoolOption
    (
        "no-mesh", // noMesh
        "Suppress writing the geometry."
        " Can be useful for converting partial results for a static geometry"
    );
    argList::addOptionCompat("no-mesh", {"noMesh", 1806});

    argList::addOption
    (
        "fields",
        "wordRes",
        "Specify single or multiple fields to write (all by default)\n"
        "Eg, 'T' or '( \"U.*\" )'"
    );

    argList::addOption
    (
        "name",
        "subdir",
        "Sub-directory name for Ensight output (default: 'Ensight')"
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

    caseOpts.width(args.get<label>("width", 8));
    caseOpts.overwrite(false); // leave existing output directory

    // Can also have separate directory for lagrangian
    // caseOpts.separateCloud(true);

    // Define sub-directory name to use for EnSight data.
    // The path to the ensight directory is at case level only
    // - For parallel cases, data only written from master
    fileName ensightDir = args.get<word>("name", "Ensight");
    if (!ensightDir.isAbsolute())
    {
        ensightDir = args.globalPath()/ensightDir;
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
    // Output configuration
    //

    // Control for renumbering iterations
    label indexingNumber = 0;
    const bool optIndex = args.readIfPresent("index", indexingNumber);
    const bool doLagrangian = !args.found("no-lagrangian");

    // Write the geometry, unless otherwise specified
    bool doGeometry = !args.found("no-mesh");

    //
    // Output configuration (field related)
    //

    wordRes fieldPatterns;
    args.readListIfPresent<wordRe>("fields", fieldPatterns);


    // Construct the list of ensight parts for the entire mesh
    ensightParts ensParts(mesh);

    // Write summary information
    if (Pstream::master())
    {
        Info<< "Converting " << timeDirs.size() << " time steps" << endl;

        OFstream info(ensCase.path()/"partsInfo");

        info
            << "// summary of ensight parts" << nl << nl;
        ensParts.writeSummary(info);
    }

    #include "checkMeshMoving.H"
    #include "findCloudFields.H"

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
        checkData(mesh, timeDirs, objectNames);

        testedObjectNames = objectNames;
    }

    if (meshMoving && !doGeometry)
    {
        Info<< "mesh is moving: ignoring '-no-mesh' option" << endl;
        doGeometry = true;
    }


    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        #include "getTimeIndex.H"
        #include "moveMesh.H"

        ensCase.setTime(timeDirs[timeI], timeIndex);

        Info<< "Time [" << timeIndex << "] = " << runTime.timeName() << nl;

        if (timeI == 0 || mesh.moving())
        {
            if (mesh.moving())
            {
                ensParts.recalculate(mesh);
            }

            if (doGeometry)
            {
                autoPtr<ensightGeoFile> os = ensCase.newGeometry(meshMoving);
                ensParts.write(os.ref());
            }
        }

        // Objects at this time
        IOobjectList objects(mesh, runTime.timeName());

        // Restrict to objects that are available for all times
        objects.filterObjects(testedObjectNames);

        // Volume, internal fields
        #include "convertVolumeFields.H"

        // Lagrangian fields
        #include "convertLagrangian.H"

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
