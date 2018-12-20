/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2018 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "IOmanip.H"
#include "Fstream.H"
#include "OSspecific.H"
#include "ensightCase.H"
#include "ensightPartFaces.H"
#include "ensightSerialOutput.H"
#include "ensightPTraits.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::ensightSurfaceWriter::writeUncollated
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const meshedSurf& surf,
    const word& fieldName,
    const Field<Type>& values,
    const bool isNodeValues,
    const bool verbose
) const
{
    const ensight::FileName surfName(surfaceName);
    const ensight::VarName  varName(fieldName);

    // geometry:  rootdir/time/<field>/surfaceName.case
    // geometry:  rootdir/time/<field>/surfaceName.<index>.mesh
    // field:     rootdir/time/<field>/surfaceName.<index>.field

    // Variable name as sub-directory for results
    // eg, something like this:
    // - VAR1/SURF1.case
    // - VAR1/SURF1.0000.mesh
    // - VAR1/SURF1.0001.VAR1
    // - VAR1/SURF2.case
    // - VAR1/SURF2.0000.mesh
    // - VAR1/SURF2.0001.VAR1
    // and
    // - VAR2/SURF1.case
    // - VAR2/SURF1.0000.mesh
    // - VAR2/SURF1.0001.VAR2

    const fileName baseDir = outputDir/varName;
    const fileName timeDir = outputDir.name();
    // Convert timeDir to a value (if possible - use 0.0 otherwise)
    scalar timeValue = 0.0;
    readScalar(timeDir, timeValue);


    if (!isDir(baseDir))
    {
        mkDir(baseDir);
    }

    OFstream osCase(baseDir/surfName + ".case", IOstream::ASCII);

    // Format options
    osCase.setf(ios_base::left);
    osCase.setf(ios_base::scientific, ios_base::floatfield);
    osCase.precision(5);

    ensightGeoFile osGeom
    (
        baseDir,
        surfName + ".00000000.mesh",
        writeFormat_
    );
    ensightFile osField
    (
        baseDir,
        surfName + ".00000000." + varName,
        writeFormat_
    );

    if (verbose)
    {
        Info<< "Writing case file to " << osCase.name() << endl;
    }

    osCase
        << "FORMAT" << nl
        << "type: ensight gold" << nl
        << nl
        << "GEOMETRY" << nl
        << "model:  1   " << osGeom.name().name() << nl
        << nl
        << "VARIABLE" << nl
        << ensightPTraits<Type>::typeName
        <<
        (
            isNodeValues
          ? " per node:    1  "  // time-set 1
          : " per element: 1  "  // time-set 1
        )
        << setw(15) << varName << ' '
        << surfName.c_str() << ".********." << varName << nl;

    osCase
        << nl
        << "TIME" << nl;

    printTimeset(osCase, 1, timeValue);
    osCase << "# end" << nl;


    ensightPartFaces ensPart
    (
        0,
        osGeom.name().name(),
        surf.points(),
        surf.faces(),
        true // contiguous points
    );
    osGeom << ensPart;

    // Write field
    osField.writeKeyword(ensightPTraits<Type>::typeName);
    ensightSerialOutput::writeField
    (
        values,
        ensPart,
        osField,
        isNodeValues
    );

    return osCase.name();
}


template<class Type>
Foam::fileName Foam::ensightSurfaceWriter::writeCollated
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const meshedSurf& surf,
    const word& fieldName,
    const Field<Type>& values,
    const bool isNodeValues,
    const bool verbose
) const
{
    const ensight::FileName surfName(surfaceName);
    const ensight::VarName  varName(fieldName);

    // geometry:  rootdir/surfaceName/surfaceName.case
    // geometry:  rootdir/surfaceName/surfaceName/data/<index>/geometry
    // field:     rootdir/surfaceName/surfaceName/data/<index>/field

    // Use surface name as sub-directory for results
    // eg, something like this:
    // - SURF1/SURF1.case
    // - SURF1/SURF1/data/00000000/geometry
    // - SURF1/SURF1/data/00000000/VAR1
    // - SURF1/SURF1/data/00000000/VAR2
    // and
    // - SURF2/SURF2.case
    // - SURF2/SURF2/data/00000000/geometry
    // - SURF2/SURF2/data/00000000/VAR1
    // - SURF2/SURF2/data/00000000/VAR2

    // Names "data" and "geometry" as per ensightCase:
    const char* fmt  = "%08d";
    const char* mask = "data/********/";

    const fileName baseDir = outputDir.path()/surfName;
    const fileName timeDir = outputDir.name();
    // Convert timeDir to a value (if possible - use 0.0 otherwise)
    scalar timeValue = 0.0;
    readScalar(timeDir, timeValue);

    scalar meshValue = 0;

    if (!isDir(baseDir))
    {
        mkDir(baseDir);
    }

    label meshIndex = 0;
    label timeIndex = 0;

    fileName geometryName;

    // Do case file
    {
        dictionary dict;
        scalarList meshes;
        scalarList times;
        bool stateChanged = false;

        if (isFile(baseDir/"fieldsDict"))
        {
            IFstream is(baseDir/"fieldsDict");
            if (is.good() && dict.read(is))
            {
                dict.readIfPresent("meshes", meshes);
                dict.readIfPresent("times", times);

                timeIndex = 1+findLower(times, timeValue);

                if (dict.readIfPresent("updateMesh", meshValue))
                {
                    meshIndex = 1+findLower(meshes, meshValue);

                    dict.remove("updateMesh");
                    stateChanged = true;
                }
                else if (meshes.size())
                {
                    meshIndex = meshes.size()-1;
                    meshValue = meshes.last();
                }
                else
                {
                    meshIndex = 0;
                }
            }
        }

        // Update stored times list
        meshes.resize(meshIndex+1, -1);
        times.resize(timeIndex+1, -1);

        if (meshes[meshIndex] != meshValue)
        {
            stateChanged = true;
        }
        if (times[timeIndex] != timeValue)
        {
            stateChanged = true;
        }

        meshes[meshIndex] = meshValue;
        times[timeIndex] = timeValue;

        geometryName =
            "data"/word::printf(fmt, meshIndex)/ensightCase::geometryName;


        // Add my information to dictionary
        {
            dict.set("meshes", meshes);
            dict.set("times", times);
            if (dict.found("fields"))
            {
                dictionary& fieldsDict = dict.subDict("fields");
                if (!fieldsDict.found(fieldName))
                {
                    dictionary fieldDict;
                    fieldDict.set("type", ensightPTraits<Type>::typeName);
                    fieldDict.set("name", varName); // ensight variable name

                    fieldsDict.set(fieldName, fieldDict);

                    stateChanged = true;
                }
            }
            else
            {
                dictionary fieldDict;
                fieldDict.set("type", ensightPTraits<Type>::typeName);
                fieldDict.set("name", varName); // ensight variable name

                dictionary fieldsDict;
                fieldsDict.set(fieldName, fieldDict);

                dict.set("fields", fieldsDict);

                stateChanged = true;
            }
        }


        if (stateChanged)
        {
            if (verbose)
            {
                Info<< "Writing state file to fieldsDict" << endl;
            }
            {
                OFstream os(baseDir/"fieldsDict");
                os << "// Summary of Ensight fields, times" << nl << nl;
                dict.write(os, false);
            }

            OFstream osCase(baseDir/surfName + ".case", IOstream::ASCII);

            // Format options
            osCase.setf(ios_base::left);
            osCase.setf(ios_base::scientific, ios_base::floatfield);
            osCase.precision(5);

            if (verbose)
            {
                Info<< "Writing case file to " << osCase.name() << endl;
            }

            // The geometry can be any of the following:
            // 0: constant/static
            // 1: moving, with the same frequency as the data
            // 2: moving, with different frequency as the data

            const label tsGeom =
                (meshes.size() == 1 ? 0 : meshes == times ? 1 : 2);

            osCase
                << "FORMAT" << nl
                << "type: ensight gold" << nl
                << nl
                << "GEOMETRY" << nl;


            if (tsGeom)
            {
                // moving
                osCase
                    << "model:  " << tsGeom << "   " // time-set (1|2)
                    << mask << geometryName.name() << nl;
            }
            else
            {
                // steady
                osCase
                    << "model:  "
                    << geometryName.c_str() << nl;
            }

            osCase
                << nl
                << "VARIABLE" << nl;

            const dictionary& fieldsDict = dict.subDict("fields");
            for (const entry& dEntry : fieldsDict)
            {
                const dictionary& subDict = dEntry.dict();

                const word fieldType(subDict.get<word>("type"));
                const word varName = subDict.lookupOrDefault
                (
                    "name",
                    dEntry.keyword() // fieldName as fallback
                );

                osCase
                    << fieldType
                    <<
                    (
                        isNodeValues
                      ? " per node:    1  "  // time-set 1
                      : " per element: 1  "  // time-set 1
                    )
                    << setw(15) << varName << ' '
                    << mask << varName << nl;
            }

            osCase
                << nl
                << "TIME" << nl;

            printTimeset(osCase, 1, times);
            if (tsGeom == 2)
            {
                printTimeset(osCase, 2, meshes);
            }

            osCase << "# end" << nl;
        }
    }


    // Location for data (and possibly the geometry as well)
    fileName dataDir = baseDir/"data"/word::printf(fmt, timeIndex);

    // As per mkdir -p "data/00000000"
    mkDir(dataDir);


    const fileName meshFile(baseDir/geometryName);

    // Write geometry
    ensightPartFaces ensPart
    (
        0,
        meshFile.name(),
        surf.points(),
        surf.faces(),
        true // contiguous points
    );
    if (!exists(meshFile))
    {
        if (verbose)
        {
            Info<< "Writing mesh file to " << meshFile.name() << endl;
        }
        // Use two-argument form for path-name to avoid validating the base-dir
        ensightGeoFile osGeom(meshFile.path(), meshFile.name(), writeFormat_);
        osGeom << ensPart;
    }

    // Write field
    ensightFile osField
    (
        dataDir,
        varName,
        writeFormat_
    );

    if (verbose)
    {
        Info<< "Writing field file to " << osField.name() << endl;
    }

    // Write field
    osField.writeKeyword(ensightPTraits<Type>::typeName);
    ensightSerialOutput::writeField
    (
        values,
        ensPart,
        osField,
        isNodeValues
    );

    // Place a timestamp in the directory for future reference
    {
        OFstream timeStamp(dataDir/"time");
        timeStamp
            << "#   timestep time" << nl
            << dataDir.name() << " " << timeValue << nl;
    }

    return baseDir/surfName + ".case";
}


template<class Type>
Foam::fileName Foam::ensightSurfaceWriter::writeTemplate
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const meshedSurf& surf,
    const word& fieldName,
    const Field<Type>& values,
    const bool isNodeValues,
    const bool verbose
) const
{
    if (collateTimes_)
    {
        return writeCollated
        (
            outputDir,
            surfaceName,
            surf,
            fieldName,
            values,
            isNodeValues,
            verbose
        );
    }
    else
    {
        return writeUncollated
        (
            outputDir,
            surfaceName,
            surf,
            fieldName,
            values,
            isNodeValues,
            verbose
        );
    }
}


// ************************************************************************* //
