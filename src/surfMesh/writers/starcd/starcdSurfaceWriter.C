/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011, 2015-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011 OpenFOAM Foundation
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

#include "starcdSurfaceWriter.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "MeshedSurfaceProxy.H"
#include "surfaceWriterMethods.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{
    defineTypeName(starcdWriter);
    addToRunTimeSelectionTable(surfaceWriter, starcdWriter, word);
}
}

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    // Emit each component
    template<class Type>
    static inline void writeData(Ostream& os, const Type& val)
    {
        const direction ncmpt = pTraits<Type>::nComponents;
        for (direction cmpt=0; cmpt < ncmpt; ++cmpt)
        {
            os  << ' ' << component(val, cmpt);
        }
        os  << nl;
    }

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriters::starcdWriter::starcdWriter()
:
    surfaceWriter()
{}


Foam::surfaceWriters::starcdWriter::starcdWriter
(
    const dictionary& options
)
:
    surfaceWriter(options)
{}


Foam::surfaceWriters::starcdWriter::starcdWriter
(
    const meshedSurf& surf,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    starcdWriter(options)
{
    open(surf, outputPath, parallel);
}


Foam::surfaceWriters::starcdWriter::starcdWriter
(
    const pointField& points,
    const faceList& faces,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    starcdWriter(options)
{
    open(points, faces, outputPath, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::surfaceWriters::starcdWriter::write()
{
    checkOpen();

    // Geometry:  rootdir/<TIME>/surfaceName.{inp,cel,vrt}

    fileName outputFile = outputPath_;
    if (useTimeDir() && !timeName().empty())
    {
        // Splice in time-directory
        outputFile = outputPath_.path() / timeName() / outputPath_.name();
    }
    outputFile.ext("inp");

    if (verbose_)
    {
        Info<< "Writing geometry to " << outputFile << endl;
    }

    const meshedSurf& surf = surface();

    if (Pstream::master() || !parallel_)
    {
        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        MeshedSurfaceProxy<face>
        (
            surf.points(),
            surf.faces()
        ).write(outputFile, "inp");
    }

    wroteGeom_ = true;
    return outputFile;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::surfaceWriters::starcdWriter::writeTemplate
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

    // Field:  rootdir/<TIME>/<field>_surfaceName.usr

    fileName outputFile = outputPath_.path();
    if (useTimeDir() && !timeName().empty())
    {
        // Splice in time-directory
        outputFile /= timeName();
    }

    // Append <field>_surfaceName.usr
    outputFile /= fieldName + '_' + outputPath_.name();
    outputFile.ext("usr");

    if (verbose_)
    {
        Info<< "Writing field " << fieldName << " to " << outputFile << endl;
    }


    // geometry merge() implicit
    tmp<Field<Type>> tfield = mergeField(localValues);

    if (Pstream::master() || !parallel_)
    {
        const auto& values = tfield();

        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        OFstream os(outputFile);

        // 1-based ids
        label elemId = 1;

        // No header, just write values
        for (const Type& val : values)
        {
            os  << elemId;
            writeData(os, val);

            ++elemId;
        }
    }

    wroteGeom_ = true;
    return outputFile;
}


// Field writing methods
defineSurfaceWriterWriteFields(Foam::surfaceWriters::starcdWriter);


// ************************************************************************* //
