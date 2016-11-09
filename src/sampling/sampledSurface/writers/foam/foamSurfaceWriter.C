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

#include "foamSurfaceWriter.H"
#include "OFstream.H"
#include "makeSurfaceWriterMethods.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceWriterType(foamSurfaceWriter);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::foamSurfaceWriter::foamSurfaceWriter()
:
    surfaceWriter()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::foamSurfaceWriter::~foamSurfaceWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::foamSurfaceWriter::write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const meshedSurf& surf,
    const bool verbose
) const
{
    const pointField& points = surf.points();
    const faceList&    faces = surf.faces();

    fileName surfaceDir(outputDir/surfaceName);

    if (!isDir(surfaceDir))
    {
        mkDir(surfaceDir);
    }

    if (verbose)
    {
        Info<< "Writing geometry to " << surfaceDir << endl;
    }


    // Points
    OFstream(surfaceDir/"points")() << points;

    // Faces
    OFstream(surfaceDir/"faces")() << faces;

    // Face centers. Not really necessary but very handy when reusing as inputs
    // for e.g. timeVaryingMapped bc.
    pointField faceCentres(faces.size(), Zero);

    forAll(faces, facei)
    {
        faceCentres[facei] = faces[facei].centre(points);
    }

    OFstream(surfaceDir/"faceCentres")() << faceCentres;

    return surfaceDir;
}


// Create write methods
defineSurfaceWriterWriteFields(Foam::foamSurfaceWriter);


// ************************************************************************* //
