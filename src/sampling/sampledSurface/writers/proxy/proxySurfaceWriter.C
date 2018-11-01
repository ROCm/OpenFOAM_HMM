/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "proxySurfaceWriter.H"
#include "MeshedSurfaceProxy.H"
#include "OSspecific.H"
#include "makeSurfaceWriterMethods.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(proxySurfaceWriter, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::proxySurfaceWriter::proxySurfaceWriter(const word& fileExt)
:
    surfaceWriter(),
    fileExtension_(fileExt)
{}


Foam::proxySurfaceWriter::proxySurfaceWriter
(
    const word& fileExt,
    const dictionary& options
)
:
    surfaceWriter(),
    fileExtension_(fileExt),
    options_(options)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::proxySurfaceWriter::write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const meshedSurf& surf,
    const bool verbose
) const
{
    // Avoid bad values
    if (fileExtension_.empty())
    {
        return fileName::null;
    }

    fileName outputFile(outputDir/surfaceName + '.' + fileExtension_);

    if (!isDir(outputFile.path()))
    {
        mkDir(outputFile.path());
    }

    if (verbose)
    {
        Info<< "Writing geometry to " << outputFile << endl;
    }

    MeshedSurfaceProxy<face>
    (
        surf.points(),
        surf.faces()
    ).write
    (
        outputFile,
        fileExtension_,
        options_
    );

    return outputFile;
}


// ************************************************************************* //
