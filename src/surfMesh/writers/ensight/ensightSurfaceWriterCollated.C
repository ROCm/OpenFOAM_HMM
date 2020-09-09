/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

int Foam::surfaceWriters::ensightWriter::geometryTimeset() const
{
    const scalarList& times = caching_.times();
    const bitSet& geoms = caching_.geometries();

    if (geoms.count() <= 1)
    {
        // Static
        return 0;
    }

    if (geoms.size() == times.size() && geoms.all())
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
    // CaseFile:  rootdir/surfaceName/surfaceName.case
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
    // Geometry changed since last output? Capture now before any merging.
    const bool geomChanged = (!upToDate_);

    checkOpen();

    const ensight::FileName surfName(outputPath_.name());
    const ensight::VarName  varName(fieldName);


    // Collated
    // ========
    // CaseFile:  rootdir/surfaceName/surfaceName.case
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



    // Implicit geometry merge()
    tmp<Field<Type>> tfield = mergeField(localValues);

    const meshedSurf& surf = surface();

    if (Pstream::master() || !parallel_)
    {
        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        const bool stateChanged =
            caching_.update
            (
                baseDir,
                timeValue,
                geomChanged,
                fieldName,
                ensightPTraits<Type>::typeName,
                varName
            );


        // The most current time and geometry indices
        const label timeIndex = caching_.latestTimeIndex();
        const label geomIndex = caching_.latestGeomIndex();


        // This will be used for the name of a static geometry,
        // or just the masking part for moving geometries.
        const fileName geometryName
        (
            "data"/word::printf(fmt, geomIndex)/ensightCase::geometryName
        );


        // Do case file
        if (stateChanged)
        {
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


            for (const entry& dEntry : caching_.fieldsDict())
            {
                const dictionary& subDict = dEntry.dict();

                const word varType(subDict.get<word>("type"));
                const word varName
                (
                    subDict.getOrDefault<word>
                    (
                        "name",
                        dEntry.keyword() // fieldName as fallback
                    )
                );

                osCase
                    << varType
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

            printTimeset(osCase, 1, caching_.times());
            if (tsGeom == 2)
            {
                printTimeset
                (
                    osCase,
                    tsGeom,
                    caching_.times(),
                    caching_.geometries()
                );
            }

            osCase << "# end" << nl;
        }


        // Location for data (and possibly the geometry as well)
        fileName dataDir = baseDir/"data"/word::printf(fmt, timeIndex);

        // As per mkdir -p "data/00000000"
        mkDir(dataDir);


        const fileName geomFile(baseDir/geometryName);

        // Ensight Geometry
        ensightOutputSurface part
        (
            surf.points(),
            surf.faces(),
            geomFile.name()
        );

        if (!exists(geomFile))
        {
            if (verbose_)
            {
                Info<< "Writing geometry to " << geomFile.name() << endl;
            }

            // Two-argument form for path-name to avoid validating base-dir
            ensightGeoFile osGeom
            (
                geomFile.path(),
                geomFile.name(),
                writeFormat_
            );
            part.write(osGeom); // serial
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

        // Write field (serial only)
        osField.writeKeyword(ensightPTraits<Type>::typeName);
        part.writeData(osField, tfield(), this->isPointData());

        // Timestamp in the directory for future reference
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
