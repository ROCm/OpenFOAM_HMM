/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "rawSurfaceWriter.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "makeSurfaceWriterMethods.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceWriterType(rawSurfaceWriter);
    addToRunTimeSelectionTable(surfaceWriter, rawSurfaceWriter, wordDict);
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    // Emit x,y,z
    static inline void writePoint(Ostream& os, const point& p)
    {
        os << p.x() << ' ' << p.y() << ' ' << p.z();
    }
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Field writing implementation
#include "rawSurfaceWriterImpl.C"

// Field writing methods
defineSurfaceWriterWriteFields(Foam::rawSurfaceWriter);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rawSurfaceWriter::rawSurfaceWriter()
:
    surfaceWriter(),
    writeCompression_(IOstream::UNCOMPRESSED)
{}


Foam::rawSurfaceWriter::rawSurfaceWriter(const dictionary& options)
:
    surfaceWriter(),
    writeCompression_(IOstream::UNCOMPRESSED)
{
    if (options.found("compression"))
    {
        writeCompression_ =
            IOstream::compressionEnum(options.get<word>("compression"));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::rawSurfaceWriter::write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const meshedSurf& surf,
    const bool verbose
) const
{
    // geometry:  rootdir/time/surfaceName.raw

    const pointField& points = surf.points();
    const faceList&    faces = surf.faces();


    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    OFstream os
    (
        outputDir/surfaceName + ".raw",
        IOstream::ASCII,
        IOstream::currentVersion,
        writeCompression_
    );

    if (verbose)
    {
        Info<< "Writing geometry to " << os.name() << endl;
    }

    // Header
    os  << "# geometry NO_DATA " << faces.size() << nl
        << "#  x  y  z" << nl;

    // Write faces centres
    for (const face& f : faces)
    {
        writePoint(os, f.centre(points));
        os << nl;
    }

    os  << nl;

    return os.name();
}


// ************************************************************************* //
