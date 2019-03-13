/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011, 2015-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "rawSurfaceWriter.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "surfaceWriterMethods.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{
    defineTypeName(rawWriter);
    addToRunTimeSelectionTable(surfaceWriter, rawWriter, word);
    addToRunTimeSelectionTable(surfaceWriter, rawWriter, wordDict);
}
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    // Emit x,y,z
    static inline void writePoint(Foam::Ostream& os, const Foam::point& p)
    {
        os << p.x() << ' ' << p.y() << ' ' << p.z();
    }

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Field writing implementation
#include "rawSurfaceWriterImpl.C"

// Field writing methods
defineSurfaceWriterWriteFields(Foam::surfaceWriters::rawWriter);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriters::rawWriter::rawWriter()
:
    surfaceWriter(),
    writeCompression_(IOstream::UNCOMPRESSED)
{}


Foam::surfaceWriters::rawWriter::rawWriter
(
    const dictionary& options
)
:
    surfaceWriter(options),
    writeCompression_
    (
        IOstream::compressionEnum
        (
            options.lookupOrDefault<word>("compression", "false")
        )
    )
{}


Foam::surfaceWriters::rawWriter::rawWriter
(
    const meshedSurf& surf,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    rawWriter(options)
{
    open(surf, outputPath, parallel);
}


Foam::surfaceWriters::rawWriter::rawWriter
(
    const pointField& points,
    const faceList& faces,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    rawWriter(options)
{
    open(points, faces, outputPath, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::surfaceWriters::rawWriter::write()
{
    checkOpen();

    // Geometry:  rootdir/<TIME>/surfaceName.raw

    fileName outputFile = outputPath_;
    if (useTimeDir() && !timeName().empty())
    {
        // Splice in time-directory
        outputFile = outputPath_.path() / timeName() / outputPath_.name();
    }
    outputFile.ext("raw");

    if (verbose_)
    {
        Info<< "Writing geometry to " << outputFile << endl;
    }


    const meshedSurf& surf = surface();

    if (Pstream::master() || !parallel_)
    {
        const pointField& points = surf.points();
        const faceList& faces = surf.faces();

        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        OFstream os
        (
            outputFile,
            IOstream::ASCII,
            IOstream::currentVersion,
            writeCompression_
        );

        // Header
        {
            os  << "# geometry NO_DATA " << faces.size() << nl
                << "#  x  y  z" << nl;
        }

        // Write faces centres
        for (const face& f : faces)
        {
            writePoint(os, f.centre(points));
            os << nl;
        }

        os  << nl;
    }

    wroteGeom_ = true;
    return outputFile;
}


// ************************************************************************* //
