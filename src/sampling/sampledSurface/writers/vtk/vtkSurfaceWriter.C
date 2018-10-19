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

#include "OFstream.H"
#include "OSspecific.H"
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
    Ostream& os,
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

    // header
    os
        << "# vtk DataFile Version 2.0" << nl
        << title.c_str() << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl;

    // Write vertex coords
    os  << "POINTS " << points.size() << " double" << nl;
    for (const point& pt : points)
    {
        os  << float(pt.x()) << ' '
            << float(pt.y()) << ' '
            << float(pt.z()) << nl;
    }
    os  << nl;


    // Write faces
    label nNodes = 0;
    for (const face& f : faces)
    {
        nNodes += f.size();
    }

    os  << "POLYGONS " << faces.size() << ' '
        << faces.size() + nNodes << nl;

    for (const face& f : faces)
    {
        os  << f.size();
        for (const label verti : f)
        {
            os  << ' ' << verti;
        }
        os  << nl;
    }
}


namespace Foam
{

    template<>
    void Foam::vtkSurfaceWriter::writeData
    (
        Ostream& os,
        const Field<scalar>& values
    )
    {
        os  << "1 " << values.size() << " double" << nl;

        forAll(values, elemI)
        {
            if (elemI)
            {
                if (elemI % 10)
                {
                    os  << ' ';
                }
                else
                {
                    os  << nl;
                }
            }

            os  << float(values[elemI]);
        }
        os  << nl;
    }


    template<>
    void Foam::vtkSurfaceWriter::writeData
    (
        Ostream& os,
        const Field<vector>& values
    )
    {
        os  << "3 " << values.size() << " double" << nl;

        for (const vector& v : values)
        {
            os  << float(v[0]) << ' '
                << float(v[1]) << ' '
                << float(v[2]) << nl;
        }
    }


    template<>
    void Foam::vtkSurfaceWriter::writeData
    (
        Ostream& os,
        const Field<sphericalTensor>& values
    )
    {
        os  << "1 " << values.size() << " double" << nl;

        for (const sphericalTensor& v : values)
        {
            os  << float(v[0]) << nl;
        }
    }


    template<>
    void Foam::vtkSurfaceWriter::writeData
    (
        Ostream& os,
        const Field<symmTensor>& values
    )
    {
        os  << "6 " << values.size() << " double" << nl;

        // symmTensor ( XX, XY, XZ, YY, YZ, ZZ )
        // VTK order  ( XX, YY, ZZ, XY, YZ, XZ ) -> (0, 3, 5, 1, 4, 2)

        for (const symmTensor& v : values)
        {
            os  << float(v[0]) << ' ' << float(v[3]) << ' ' << float(v[5])
                << ' '
                << float(v[1]) << ' ' << float(v[4]) << ' ' << float(v[2])
                << nl;
        }
    }


    template<>
    void Foam::vtkSurfaceWriter::writeData
    (
        Ostream& os,
        const Field<tensor>& values
    )
    {
        os  << "9 " << values.size() << " double" << nl;

        for (const tensor& v : values)
        {
            os  << float(v[0]) << ' ' << float(v[1]) << ' ' << float(v[2])
                << ' '
                << float(v[3]) << ' ' << float(v[4]) << ' ' << float(v[5])
                << ' '
                << float(v[6]) << ' ' << float(v[7]) << ' ' << float(v[8])
                << nl;
        }
    }
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

    opts.ascii
    (
        options.found("format")
     && (IOstream::formatEnum(options.get<word>("format")) == IOstream::ASCII)
    );

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

    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    OFstream os(outputDir/surfaceName + ".vtk");
    os.precision(precision_);

    if (verbose)
    {
        Info<< "Writing geometry to " << os.name() << endl;
    }

    writeGeometry(os, surf, surfaceName);

    return os.name();
}


// ************************************************************************* //
