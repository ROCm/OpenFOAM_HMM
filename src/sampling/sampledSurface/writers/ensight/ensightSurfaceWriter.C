/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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
    writeFormat_(IOstream::ASCII),
    collateTimes_(true)
{
    // choose ascii or binary format
    if (options.found("format"))
    {
        writeFormat_ = IOstream::formatEnum(options.lookup("format"));
    }
    options.readIfPresent("collateTimes", collateTimes_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightSurfaceWriter::~ensightSurfaceWriter()
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


Foam::fileName Foam::ensightSurfaceWriter::write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const meshedSurf& surf,
    const bool verbose
) const
{
    const pointField& points = surf.points();
    const faceList&   faces  = surf.faces();
    const ensight::FileName surfName(surfaceName);

    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    const scalar timeValue = 0.0;

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

    osCase
        << "FORMAT" << nl
        << "type: ensight gold" << nl
        << nl
        << "GEOMETRY" << nl
        << "model:        1     " << osGeom.name().name() << nl
        << nl
        << "TIME" << nl
        << "time set:                      1" << nl
        << "number of steps:               1" << nl
        << "filename start number:         0" << nl
        << "filename increment:            1" << nl
        << "time values:" << nl
        << "    " << timeValue << nl
        << nl;

    ensightPartFaces ensPart(0, osGeom.name().name(), points, faces, true);
    osGeom << ensPart;

    return osCase.name();
}


// create write methods
defineSurfaceWriterWriteFields(Foam::ensightSurfaceWriter);


// ************************************************************************* //
