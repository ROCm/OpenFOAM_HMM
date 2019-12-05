/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2015-2019 OpenCFD Ltd.
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
// The returned index is always 0 or larger (no negative values).
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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::surfaceWriters::ensightWriter::readPreviousTimes
(
    const fileName& baseDir,
    const word& dictName,
    const scalar& timeValue
)
{
    // In 1906 and earlier, the fieldsDict contained "meshes" and "times"
    // entries, each with their own time values.
    // This makes it more difficult to define the exact correspondence
    // between geometry intervals and times.
    //
    // We now instead track used geometry intervals as a bitSet.


    // Only called from master

    label timeIndex = 0;

    labelList geomIndices;
    scalarList meshTimes;

    cache_.clear();

    const fileName dictFile(baseDir/dictName);

    if (isFile(dictFile))
    {
        IFstream is(dictFile);

        if (is.good() && cache_.read(is))
        {
            meshes_.clear();

            cache_.readIfPresent("times", times_);
            timeIndex = findTimeIndex(times_, timeValue);

            if (cache_.readIfPresent("geometry", geomIndices))
            {
                // Convert indices to bitSet entries
                meshes_.set(geomIndices);
            }
            else if (cache_.readIfPresent("meshes", meshTimes))
            {
                WarningInFunction
                    << nl
                    << "Setting geometry timeset information from time values"
                    << " (fieldsDict from an older OpenFOAM version)." << nl
                    << "This may not be fully reliable." << nl
                    << nl;

                for (const scalar& meshTime : meshTimes)
                {
                    const label meshIndex = findTimeIndex(times_, meshTime);
                    meshes_.set(meshIndex);
                }
            }

            // Make length consistent with time information.
            // We read/write the indices instead of simply dumping the bitSet.
            // This makes the contents more human readable.
            meshes_.resize(times_.size());
        }
    }

    return timeIndex;
}


int Foam::surfaceWriters::ensightWriter::geometryTimeset() const
{
    if (meshes_.count() <= 1)
    {
        // Static
        return 0;
    }

    if (meshes_.size() == times_.size() && meshes_.all())
    {
        // Geometry changing is the same as fields changing
        return 1;
    }

    // Geometry changing differently from fields
    return 2;
}


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

        bool stateChanged = meshChanged;

        const label timeIndex =
        (
            times_.empty()
          ? readPreviousTimes(baseDir, "fieldsDict", timeValue)
          : findTimeIndex(times_, timeValue)
        );


        // Update stored times list and mesh index

        if (timeIndex < meshes_.size()-1)
        {
            // Clear old content when shrinking
            meshes_.unset(timeIndex);
        }

        // Extend or truncate list
        meshes_.resize(timeIndex+1);
        times_.resize(timeIndex+1, VGREAT);

        if (meshChanged)
        {
            meshes_.set(timeIndex);
        }

        if (!equalTimes(times_[timeIndex], timeValue))
        {
            stateChanged = true;
            times_[timeIndex] = timeValue;
        }


        // The most current geometry index
        const label geomIndex(max(0, meshes_.find_last()));


        // This will be used for the name of a static geometry,
        // or just the masking part for moving geometries.
        const fileName geometryName
        (
            "data"/word::printf(fmt, geomIndex)/ensightCase::geometryName
        );


        // Do case file
        {
            // Add time information to dictionary
            cache_.set("geometry", meshes_.sortedToc());
            cache_.set("times", times_);

            // Debugging, or if needed for older versions:
            //// cache_.set
            //// (
            ////     "meshes",
            ////     IndirectList<scalar>(times_, meshes_.sortedToc())
            //// );

            // Add field information to dictionary
            dictionary& fieldsDict = cache_.subDictOrAdd("fields");
            dictionary& fieldDict = fieldsDict.subDictOrAdd(fieldName);

            if (fieldDict.empty())
            {
                fieldDict.set("type", ensightPTraits<Type>::typeName);
                fieldDict.set("name", varName); // ensight variable name
                stateChanged = true;
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
                    cache_.write(os, false);
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

                const label tsGeom = geometryTimeset();

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


                for (const entry& dEntry : fieldsDict)
                {
                    const dictionary& subDict = dEntry.dict();

                    const word fieldType(subDict.get<word>("type"));
                    const word varName
                    (
                        subDict.getOrDefault<word>
                        (
                            "name",
                            dEntry.keyword() // fieldName as fallback
                        )
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

                printTimeset(osCase, 1, times_);
                if (tsGeom == 2)
                {
                    printTimeset(osCase, 2, times_, meshes_);
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
                << dataDir.name() << ' ' << timeValue << nl;
        }
    }

    wroteGeom_ = true;
    return outputFile;
}


// ************************************************************************* //
