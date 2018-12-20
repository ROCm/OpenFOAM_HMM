/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "ensightSurfaceWriter.H"
#include "ensightPartFaces.H"
#include "makeSurfaceWriterMethods.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceWriterType(ensightSurfaceWriter);
    addToRunTimeSelectionTable(surfaceWriter, ensightSurfaceWriter, wordDict);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Field writing implementation
#include "ensightSurfaceWriterImpl.C"

// Field writing methods
defineSurfaceWriterWriteFields(Foam::ensightSurfaceWriter);


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::ensightSurfaceWriter::printTimeset
(
    OSstream& os,
    const label ts,
    const scalar& timeValue
)
{
    os
        << "time set:               " << ts << nl
        << "number of steps:        " << 1 << nl;

    // Assume to be contiguous numbering
    os  << "filename start number:  0" << nl
        << "filename increment:     1" << nl
        << "time values:" << nl;

    os  << "    " << timeValue
        << nl << nl;
}


void Foam::ensightSurfaceWriter::printTimeset
(
    OSstream& os,
    const label ts,
    const UList<scalar>& values
)
{
    label count = values.size();
    os
        << "time set:               " << ts << nl
        << "number of steps:        " << count << nl;

    // Assume to be contiguous numbering
    os  << "filename start number:  0" << nl
        << "filename increment:     1" << nl
        << "time values:" << nl;

    count = 0;
    for (const scalar& t : values)
    {
        os << ' ' << setw(12) << t;

        if (++count % 6 == 0)
        {
            os << nl;
        }
    }
    os  << nl << nl;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightSurfaceWriter::ensightSurfaceWriter()
:
    surfaceWriter(),
    writeFormat_(IOstream::ASCII),
    collateTimes_(true)
{}


Foam::ensightSurfaceWriter::ensightSurfaceWriter(const dictionary& options)
:
    surfaceWriter(),
    writeFormat_
    (
        IOstreamOption::formatNames.lookupOrDefault
        (
            "format",
            options,
            IOstreamOption::ASCII,
            true  // Failsafe behaviour
        )
    ),
    collateTimes_(options.lookupOrDefault("collateTimes", true))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Note that ensight does supports geometry in a separate file,
// but setting this true leaves mesh files in the wrong places
// (when there are fields).
//
// Make this false to let the field writers take back control
bool Foam::ensightSurfaceWriter::separateGeometry() const
{
    return false;
}


void Foam::ensightSurfaceWriter::updateMesh
(
    const fileName& outputDir,
    const fileName& surfaceName
) const
{
    if (collateTimes_ && Pstream::master())
    {
        const ensight::FileName surfName(surfaceName);

        const fileName baseDir = outputDir.path()/surfName;
        const fileName timeDir = outputDir.name();

        // Convert timeDir to a value (if possible - use 0.0 otherwise)
        scalar timeValue = 0.0;
        readScalar(timeDir, timeValue);

        if (!isDir(baseDir))
        {
            mkDir(baseDir);
        }

        dictionary dict;

        if (isFile(baseDir/"fieldsDict"))
        {
            IFstream is(baseDir/"fieldsDict");
            if (is.good() && dict.read(is))
            {
                // ... any futher actions
            }
        }

        dict.set("updateMesh", timeValue);

        {
            OFstream os(baseDir/"fieldsDict");
            os << "// Summary of Ensight fields, times" << nl << nl;
            dict.write(os, false);
        }
    }
}


Foam::fileName Foam::ensightSurfaceWriter::write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const meshedSurf& surf,
    const bool verbose
) const
{
    const ensight::FileName surfName(surfaceName);

    // Uncollated
    // ==========
    // geometry:  rootdir/time/surfaceName.case
    // geometry:  rootdir/time/surfaceName.00000000.mesh

    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    OFstream osCase(outputDir/surfName + ".case");
    ensightGeoFile osGeom
    (
        outputDir,
        surfName + ".00000000.mesh",
        writeFormat_
    );

    if (verbose)
    {
        Info<< "Writing case file to " << osCase.name() << endl;
    }

    const pointField& points = surf.points();
    const faceList&   faces  = surf.faces();

    osCase
        << "FORMAT" << nl
        << "type: ensight gold" << nl
        << nl
        << "GEOMETRY" << nl
        << "model:        1     " << osGeom.name().name() << nl
        << nl
        << "TIME" << nl;

    printTimeset(osCase, 1, 0.0);

    ensightPartFaces ensPart(0, osGeom.name().name(), points, faces, true);
    osGeom << ensPart;

    return osCase.name();


    // Collated?
    // ========
    // geometry:  rootdir/surfaceName/surfaceName.case
    // geometry:  rootdir/surfaceName/surfaceName.mesh
}


// ************************************************************************* //
