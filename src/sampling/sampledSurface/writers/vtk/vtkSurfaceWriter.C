/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "vtkSurfaceWriter.H"
#include "foamVtkOutputOptions.H"
#include "OSspecific.H"
#include <fstream>
#include "makeSurfaceWriterMethods.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceWriterType(vtkSurfaceWriter);
    addToRunTimeSelectionTable(surfaceWriter, vtkSurfaceWriter, wordDict);
}

// Field writing implementation
#include "vtkSurfaceWriterImpl.C"

// Field writing methods
defineSurfaceWriterWriteFields(Foam::vtkSurfaceWriter);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtkSurfaceWriter::writeGeometry
(
    vtk::formatter& format,
    const meshedSurf& surf,
    std::string title
)
{
    const pointField& points = surf.points();
    const faceList&    faces = surf.faces();

    if (title.empty())
    {
        title = "sampleSurface";
    }

    vtk::legacy::fileHeader
    (
        format,
        title,
        vtk::fileTag::POLY_DATA
    );

    vtk::legacy::beginPoints(format.os(), points.size());

    vtk::writeList(format, points);
    format.flush();

    // Write faces
    // connectivity count without additional storage (done internally)
    label nConnectivity = 0;
    for (const face& f : faces)
    {
        nConnectivity += f.size();
    }

    vtk::legacy::beginPolys
    (
        format.os(),
        faces.size(),
        nConnectivity
    );

    for (const face& f : faces)
    {
        format.write(f.size());  // The size prefix
        vtk::writeList(format, f);
    }

    format.flush();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtkSurfaceWriter::vtkSurfaceWriter()
:
    surfaceWriter(),
    fmtType_(unsigned(vtk::formatType::LEGACY_ASCII)),
    precision_(IOstream::defaultPrecision())
{}


Foam::vtkSurfaceWriter::vtkSurfaceWriter(const dictionary& options)
:
    surfaceWriter(),
    fmtType_(static_cast<unsigned>(vtk::formatType::LEGACY_ASCII)),
    precision_(IOstream::defaultPrecision())
{
#if 0
    // Future
    // format: ascii | binary
    // legacy  true | false

    vtk::outputOptions opts(static_cast<vtk::formatType>(fmtType_));

    const word formatName = options.lookupOrDefault<word>("format", "");
    if (formatName.size())
    {
        opts.ascii
        (
            IOstream::formatEnum(formatName) == IOstream::ASCII
        );
    }

    if (options.lookupOrDefault("legacy", false))
    {
        opts.legacy(true);
    }

    // Convert back to raw data type
    fmtType_ = static_cast<unsigned>(opts.fmt());
#endif

    // The write precision for ASCII formatters
    precision_ =
        options.lookupOrDefaultCompat
        (
            "precision", {{"writePrecision", -1806}},
            IOstream::defaultPrecision()
        );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::vtkSurfaceWriter::write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const meshedSurf& surf,
    const bool verbose
) const
{
    // geometry:  rootdir/time/surfaceName.{vtk|vtp}

    fileName outputFile(outputDir/surfaceName + ".vtk");

    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }
    if (verbose)
    {
        Info<< "Writing geometry to " << outputFile << endl;
    }

    // As vtk::outputOptions
    vtk::outputOptions opts(static_cast<vtk::formatType>(fmtType_));
    opts.legacy(true);
    opts.precision(precision_);

    std::ofstream os(outputFile);

    auto format = opts.newFormatter(os);

    writeGeometry(*format, surf, surfaceName);

    return outputFile;
}


// ************************************************************************* //
