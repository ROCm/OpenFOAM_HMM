/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "boundaryDataSurfaceWriter.H"
#include "argList.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "IOmanip.H"
#include "Time.H"
#include "pointIOField.H"
#include "primitivePatch.H"
#include "surfaceWriterMethods.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{
    defineTypeName(boundaryDataWriter);
    addToRunTimeSelectionTable(surfaceWriter, boundaryDataWriter, word);
    addToRunTimeSelectionTable(surfaceWriter, boundaryDataWriter, wordDict);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriters::boundaryDataWriter::boundaryDataWriter()
:
    surfaceWriter(),
    header_(true),
    streamOpt_(),
    fieldScale_()
{}


Foam::surfaceWriters::boundaryDataWriter::boundaryDataWriter
(
    const dictionary& options
)
:
    surfaceWriter(options),
    header_(options.getOrDefault("header", true)),
    streamOpt_
    (
        IOstream::formatEnum("format", options, IOstream::ASCII),
        IOstream::compressionEnum("compression", options)
    ),
    fieldScale_(options.subOrEmptyDict("fieldScale"))
{}


Foam::surfaceWriters::boundaryDataWriter::boundaryDataWriter
(
    const meshedSurf& surf,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    boundaryDataWriter(options)
{
    open(surf, outputPath, parallel);
}


Foam::surfaceWriters::boundaryDataWriter::boundaryDataWriter
(
    const pointField& points,
    const faceList& faces,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    boundaryDataWriter(options)
{
    open(points, faces, outputPath, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfaceWriters::boundaryDataWriter::serialWriteGeometry
(
    const regIOobject& iopts,
    const meshedSurf& surf
)
{
    const pointField& points = surf.points();
    const faceList& faces = surf.faces();

    if (verbose_)
    {
        if (this->isPointData())
        {
            Info<< "Writing points: " << iopts.objectPath() << endl;
        }
        else
        {
            Info<< "Writing face centres: " << iopts.objectPath() << endl;
        }
    }

    // Like regIOobject::writeObject without instance() adaptation
    // since this would write to e.g. 0/ instead of postProcessing/

    OFstream osGeom(iopts.objectPath(), streamOpt_);

    if (header_)
    {
        iopts.writeHeader(osGeom);
    }

    if (this->isPointData())
    {
        // Just like writeData, but without copying beforehand
        osGeom << points;
    }
    else
    {
        primitivePatch pp(SubList<face>(faces), points);

        // Just like writeData, but without copying beforehand
        osGeom << pp.faceCentres();
    }

    if (header_)
    {
        iopts.writeEndDivider(osGeom);
    }
}


Foam::fileName Foam::surfaceWriters::boundaryDataWriter::write()
{
    checkOpen();

    // Geometry: rootdir/surfaceName/"points"
    // Field:    rootdir/surfaceName/<TIME>/field

    fileName surfaceDir = outputPath_;

    // Dummy Time to use as objectRegistry
    autoPtr<Time> dummyTimePtr(Time::New(argList::envGlobalPath()));

    const meshedSurf& surf = surface();

    if (Pstream::master() || !parallel_)
    {
        if (!isDir(surfaceDir))
        {
            mkDir(surfaceDir);
        }

        // Write sample locations
        pointIOField iopts
        (
            IOobject
            (
                surfaceDir/"points",
                *dummyTimePtr,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );
        iopts.note() = (this->isPointData() ? "point data" : "face data");

        serialWriteGeometry(iopts, surf);
    }

    wroteGeom_ = true;
    return surfaceDir;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::surfaceWriters::boundaryDataWriter::writeTemplate
(
    const word& fieldName,
    const Field<Type>& localValues
)
{
    checkOpen();

    // Geometry: rootdir/surfaceName/"points"
    // Field:    rootdir/surfaceName/<TIME>/field

    fileName surfaceDir = outputPath_;

    const fileName outputFile(surfaceDir/timeName()/fieldName);


    // Output scaling for the variable, but not for integer types.
    // could also solve with clever templating

    const scalar varScale =
    (
        std::is_integral<Type>::value
      ? scalar(1)
      : fieldScale_.getOrDefault<scalar>(fieldName, 1)
    );


    // Dummy Time to use as objectRegistry
    autoPtr<Time> dummyTimePtr(Time::New(argList::envGlobalPath()));


    // Implicit geometry merge()
    tmp<Field<Type>> tfield = mergeField(localValues) * varScale;

    const meshedSurf& surf = surface();

    if (Pstream::master() || !parallel_)
    {
        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        // Write sample locations
        {
            pointIOField iopts
            (
                IOobject
                (
                    surfaceDir/"points",
                    *dummyTimePtr,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );
            iopts.note() = (this->isPointData() ? "point data" : "face data");

            serialWriteGeometry(iopts, surf);
        }

        // Write field
        {
            IOField<Type> iofld
            (
                IOobject
                (
                    outputFile,
                    *dummyTimePtr,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );
            iofld.note() = (this->isPointData() ? "point data" : "face data");

            OFstream osField(iofld.objectPath(), streamOpt_);

            if (header_)
            {
                iofld.writeHeader(osField);
            }

            // Just like writeData, but without copying beforehand
            osField << tfield();

            if (header_)
            {
                iofld.writeEndDivider(osField);
            }
        }
    }

    wroteGeom_ = true;
    return surfaceDir;
}


// Field writing methods
defineSurfaceWriterWriteFields(Foam::surfaceWriters::boundaryDataWriter);


// ************************************************************************* //
