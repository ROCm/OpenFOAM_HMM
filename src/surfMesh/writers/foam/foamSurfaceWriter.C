/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2019 OpenCFD Ltd.
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

#include "foamSurfaceWriter.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "surfaceWriterMethods.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{
    defineTypeName(foamWriter);
    addToRunTimeSelectionTable(surfaceWriter, foamWriter, word);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriters::foamWriter::foamWriter()
:
    surfaceWriter()
{}


Foam::surfaceWriters::foamWriter::foamWriter
(
    const dictionary& options
)
:
    surfaceWriter(options)
{}


Foam::surfaceWriters::foamWriter::foamWriter
(
    const meshedSurf& surf,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    foamWriter(options)
{
    open(surf, outputPath, parallel);
}


Foam::surfaceWriters::foamWriter::foamWriter
(
    const pointField& points,
    const faceList& faces,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    foamWriter(options)
{
    open(points, faces, outputPath, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::surfaceWriters::foamWriter::write()
{
    checkOpen();

    // Geometry:
    // - rootdir/<TIME>/surfaceName/{points,faces}

    fileName surfaceDir = outputPath_;
    if (useTimeDir() && !timeName().empty())
    {
        // Splice in time-directory
        surfaceDir = outputPath_.path() / timeName() / outputPath_.name();
    }

    if (verbose_)
    {
        Info<< "Writing geometry to " << surfaceDir << endl;
    }


    const meshedSurf& surf = surface();

    if (Pstream::master() || !parallel_)
    {
        const pointField& points = surf.points();
        const faceList& faces = surf.faces();

        if (!isDir(surfaceDir))
        {
            mkDir(surfaceDir);
        }

        // Points
        OFstream(surfaceDir/"points")() << points;

        // Faces
        OFstream(surfaceDir/"faces")() << faces;

        // Face centers.
        // Not really necessary but very handy when reusing as inputs
        // for e.g. timeVaryingMapped bc.
        pointField faceCentres(faces.size(), Zero);

        forAll(faces, facei)
        {
            faceCentres[facei] = faces[facei].centre(points);
        }

        OFstream(surfaceDir/"faceCentres")() << faceCentres;
    }

    wroteGeom_ = true;
    return surfaceDir;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::surfaceWriters::foamWriter::writeTemplate
(
    const word& fieldName,
    const Field<Type>& localValues
)
{
    // Separate geometry
    if (!wroteGeom_)
    {
        write();
    }

    checkOpen();

    // Geometry should already have been written
    // Values to separate directory (e.g. "scalarField/p")

    // Field:    rootdir/<TIME>/surfaceName/fieldType/field
    //?? -> or        rootdir/surfaceName/fieldType/<TIME>/field

    fileName surfaceDir = outputPath_;
    if (useTimeDir() && !timeName().empty())
    {
        // Splice in time-directory
        surfaceDir = outputPath_.path() / timeName() / outputPath_.name();
    }

    const fileName outputFile
    (
        surfaceDir
      / (word(pTraits<Type>::typeName) + FieldBase::typeName)
      / fieldName
    );

    if (verbose_)
    {
        Info<< "Writing field " << fieldName << " to " << surfaceDir << endl;
    }


    // geometry merge() implicit
    tmp<Field<Type>> tfield =  mergeField(localValues);

    if (Pstream::master())
    {
        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        // Write field
        OFstream(outputFile)() << tfield();
    }

    wroteGeom_ = true;
    return outputFile;
}


// Field writing methods
defineSurfaceWriterWriteFields(Foam::surfaceWriters::foamWriter);


// ************************************************************************* //
