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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::fileName Foam::surfaceWriters::ensightWriter::writeUncollated()
{
    checkOpen();

    const ensight::FileName surfName(outputPath_.name());


    // Uncollated
    // ==========
    // CaseFile:  rootdir/<TIME>/surfaceName.case
    // Geometry:  rootdir/<TIME>/surfaceName.00000000.mesh

    fileName outputDir;
    if (useTimeDir() && !timeName().empty())
    {
        // Splice in time-directory
        outputDir = outputPath_.path() / timeName();
    }
    else
    {
        outputDir = outputPath_.path();
    }

    const fileName outputFile = outputDir / surfName + ".case";

    if (verbose_)
    {
        Info<< "Writing case file to " << outputFile << endl;
    }


    const meshedSurf& surf = surface();

    if (Pstream::master() || !parallel_)
    {
        if (!isDir(outputDir))
        {
            mkDir(outputDir);
        }

        OFstream osCase(outputFile);
        ensightGeoFile osGeom
        (
            outputDir,
            surfName + ".00000000.mesh",
            writeFormat_
        );

        osCase
            << "FORMAT" << nl
            << "type: ensight gold" << nl
            << nl
            << "GEOMETRY" << nl
            << "model:        1     " << osGeom.name().name() << nl
            << nl
            << "TIME" << nl;

        printTimeset(osCase, 1, scalar(0));

        ensightOutputSurface part
        (
            surf.points(),
            surf.faces(),
            osGeom.name().name()
        );
        part.write(osGeom); // serial
    }

    wroteGeom_ = true;
    return outputFile;
}


template<class Type>
Foam::fileName Foam::surfaceWriters::ensightWriter::writeUncollated
(
    const word& fieldName,
    const Field<Type>& localValues
)
{
    checkOpen();

    const ensight::FileName surfName(outputPath_.name());
    const ensight::VarName  varName(fieldName);


    // Uncollated
    // ==========
    // CaseFile:  rootdir/time/<field>/surfaceName.case
    // Geometry:  rootdir/time/<field>/surfaceName.<index>.mesh
    // Field:     rootdir/time/<field>/surfaceName.<index>.<field>

    // Variable name as sub-directory for results. Eg,
    // - VAR1/SURF1.case
    // - VAR1/SURF1.00000000.mesh
    // - VAR1/SURF1.00000001.VAR1
    // and
    // - VAR2/SURF1.case
    // - VAR2/SURF1.00000000.mesh
    // - VAR2/SURF1.00000001.VAR2

    fileName outputDir;
    if (useTimeDir() && !timeName().empty())
    {
        // Splice in time-directory
        outputDir = outputPath_.path() / timeName();
    }
    else
    {
        outputDir = outputPath_.path();
    }

    const fileName baseDir = outputDir / varName;
    const word   timeDir = timeName();
    const scalar timeValue = currTime_.value();

    const fileName outputFile = baseDir / surfName + ".case";

    if (verbose_)
    {
        Info<< "Writing case file to " << outputFile << endl;
    }


    // geometry merge() implicit
    tmp<Field<Type>> tfield = mergeField(localValues);

    const meshedSurf& surf = surface();

    if (Pstream::master() || !parallel_)
    {
        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        OFstream osCase(outputFile, IOstream::ASCII);

        // Format options
        osCase.setf(ios_base::left);
        osCase.setf(ios_base::scientific, ios_base::floatfield);
        osCase.precision(5);

        // Two-argument form for path-name to avoid validating base-dir
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
                this->isPointData()
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


        // Ensight Geometry
        ensightOutputSurface part
        (
            surf.points(),
            surf.faces(),
            osGeom.name().name()
        );
        part.write(osGeom); // serial

        // Write field (serial)
        osField.writeKeyword(ensightPTraits<Type>::typeName);
        part.writeData(osField, tfield(), this->isPointData());
    }

    wroteGeom_ = true;
    return outputFile;
}


// ************************************************************************* //
