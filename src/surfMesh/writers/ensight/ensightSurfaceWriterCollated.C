/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2014 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Compare time values with tolerance
static const equalOp<scalar> equalTimes(ROOTSMALL);

// Use ListOps findLower (with tolerance), to find the location of the next
// time-related index.
static label findTimeIndex(const UList<scalar>& list, const scalar val)
{
    label idx =
        findLower
        (
            list,
            val,
            0,
            [](const scalar a, const scalar b)
            {
                return (a < b) && (Foam::mag(b - a) > ROOTSMALL);
            }
        );

    if (idx < 0 || !equalTimes(list[idx], val))
    {
        ++idx;
    }

    return idx;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::surfaceWriters::ensightWriter::writeCollated()
{
    // Collated?
    // ========
    // Geometry:  rootdir/surfaceName/surfaceName.case
    // Geometry:  rootdir/surfaceName/surfaceName.mesh

    wroteGeom_ = true;
    return fileName::null;
}


template<class Type>
Foam::fileName Foam::surfaceWriters::ensightWriter::writeCollated
(
    const word& fieldName,
    const Field<Type>& localValues
)
{
    checkOpen();

    const ensight::FileName surfName(outputPath_.name());
    const ensight::VarName  varName(fieldName);


    // Collated
    // ========
    // Geometry:  rootdir/surfaceName/surfaceName.case
    // Geometry:  rootdir/surfaceName/data/<index>/geometry
    // Field:     rootdir/surfaceName/data/<index>/field

    // Use surface name as sub-directory for results. Eg,
    // - SURF1/SURF1.case
    // - SURF1/data/00000000/geometry
    // - SURF1/data/00000000/VAR1
    // - SURF1/data/00000000/VAR2

    // Names "data" and "geometry" as per ensightCase:
    const char* fmt  = "%08d";
    const char* mask = "data/********/";


    // Ignore the useTimeDir setting - manage ourselves
    const fileName baseDir = outputPath_;

    const word   timeDir = timeName();
    const scalar timeValue = currTime_.value();

    const fileName outputFile = baseDir / surfName + ".case";

    if (verbose_)
    {
        Info<< "Writing case file to " << outputFile << endl;
    }


    // Mesh changed since last output? Do before any merging.
    const bool meshChanged = (!upToDate_);

    // geometry merge() implicit
    tmp<Field<Type>> tfield = mergeField(localValues);

    const meshedSurf& surf = surface();

    if (Pstream::master() || !parallel_)
    {
        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        scalar meshValue = timeValue;
        label meshIndex = 0;
        label timeIndex = 0;

        fileName geometryName;

        // Do case file
        {
            dictionary dict;
            scalarList meshes;
            scalarList times;
            bool stateChanged = meshChanged;

            if (isFile(baseDir/"fieldsDict"))
            {
                IFstream is(baseDir/"fieldsDict");
                if (is.good() && dict.read(is))
                {
                    dict.readIfPresent("meshes", meshes);
                    dict.readIfPresent("times", times);

                    timeIndex = findTimeIndex(times, timeValue);

                    if (meshChanged)
                    {
                        meshValue = timeValue;
                        meshIndex = findTimeIndex(meshes, meshValue);
                    }
                    else if (meshes.size())
                    {
                        meshIndex = meshes.size()-1;
                        meshValue = meshes.last();
                    }
                }
            }

            // Update stored times list
            meshes.resize(meshIndex+1, -1);
            times.resize(timeIndex+1, -1);

            if
            (
                !equalTimes(meshes[meshIndex], meshValue)
             || !equalTimes(times[timeIndex], timeValue)
            )
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
                if (verbose_)
                {
                    Info<< "Writing state file to fieldsDict" << endl;
                }
                {
                    OFstream os(baseDir/"fieldsDict");
                    os << "// Summary of Ensight fields, times" << nl << nl;
                    dict.write(os, false);
                }

                OFstream osCase(outputFile, IOstream::ASCII);

                // Format options
                osCase.setf(ios_base::left);
                osCase.setf(ios_base::scientific, ios_base::floatfield);
                osCase.precision(5);

                if (verbose_)
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
                            this->isPointData()
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
            if (verbose_)
            {
                Info<< "Writing mesh file to " << meshFile.name() << endl;
            }
            // Use two-argument form for path-name to avoid validating base-dir
            ensightGeoFile osGeom
            (
                meshFile.path(),
                meshFile.name(),
                writeFormat_
            );
            osGeom << ensPart;
        }

        // Write field
        ensightFile osField
        (
            dataDir,
            varName,
            writeFormat_
        );

        if (verbose_)
        {
            Info<< "Writing field file to " << osField.name() << endl;
        }

        // Write field
        osField.writeKeyword(ensightPTraits<Type>::typeName);

        if (this->isPointData())
        {
            ensightOutput::Serial::writePointField
            (
                tfield(),
                ensPart,
                osField
                // serial
            );
        }
        else
        {
            ensightOutput::Detail::writeFaceField
            (
                tfield(),
                ensPart,
                osField,
                false // serial
            );
        }

        // Place a timestamp in the directory for future reference
        {
            OFstream timeStamp(dataDir/"time");
            timeStamp
                << "#   timestep time" << nl
                << dataDir.name() << " " << timeValue << nl;
        }
    }

    wroteGeom_ = true;
    return outputFile;
}


// ************************************************************************* //
