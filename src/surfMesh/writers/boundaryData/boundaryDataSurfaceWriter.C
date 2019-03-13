/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2015 OpenFOAM Foundation
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
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriters::boundaryDataWriter::boundaryDataWriter()
:
    surfaceWriter()
{}


Foam::surfaceWriters::boundaryDataWriter::boundaryDataWriter
(
    const dictionary& options
)
:
    surfaceWriter(options)
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

Foam::fileName Foam::surfaceWriters::boundaryDataWriter::write()
{
    checkOpen();

    // Geometry: rootdir/surfaceName/"points"
    // Field:    rootdir/surfaceName/<TIME>/field

    fileName surfaceDir = outputPath_;

    // Write points
    if (verbose_)
    {
        Info<< "Writing points to " << surfaceDir/"points" << endl;
    }


    // Dummy time to use as an objectRegistry
    const fileName caseDir(argList::envGlobalPath());

    Time dummyTime
    (
        caseDir.path(), // root-path,
        caseDir.name(), // case-name,
        "system",       //
        "constant",     //
        false,          // no function objects
        false           // no libs
    );


    const meshedSurf& surf = surface();

    if (Pstream::master() || !parallel_)
    {
        if (!isDir(surfaceDir))
        {
            mkDir(surfaceDir);
        }

        pointIOField pts
        (
            IOobject
            (
                surfaceDir/"points",
                dummyTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            surf.points()
        );

        // Do like regIOobject::writeObject but don't do instance() adaptation
        // since this would write to e.g. 0/ instead of postProcessing/

        // Try opening an OFstream for object
        OFstream os(pts.objectPath());

        //pts.writeHeader(os);
        pts.writeData(os);
        //pts.writeEndDivider(os);
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


    // Dummy time to use as an objectRegistry
    const fileName caseDir(argList::envGlobalPath());

    Time dummyTime
    (
        caseDir.path(), // root-path,
        caseDir.name(), // case-name,
        "system",       //
        "constant",     //
        false,          // no function objects
        false           // no libs
    );


    // Geometry merge() implicit
    tmp<Field<Type>> tfield =  mergeField(localValues);

    const meshedSurf& surf = surface();

    if (Pstream::master() || !parallel_)
    {
        const pointField& points = surf.points();
        const faceList& faces = surf.faces();

        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        pointIOField pts
        (
            IOobject
            (
                surfaceDir/"points",
                dummyTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            label(0)
        );

        if (this->isPointData())
        {
            if (verbose_)
            {
                Info<< "Writing points to "
                    << surfaceDir/"points" << endl;
            }
            pts = points;
        }
        else
        {
            if (verbose_)
            {
                Info<< "Writing face centres to "
                    << surfaceDir/"points" << endl;
            }

            primitivePatch pp(SubList<face>(faces, faces.size()), points);

            pts = pp.faceCentres();
        }

        {
            // Do like regIOobject::writeObject but don't do instance()
            // adaptation
            // since this would write to e.g. 0/ instead of postProcessing/

            // Try opening an OFstream for object
            OFstream os(pts.objectPath());

            //pts.writeHeader(os);
            pts.writeData(os);
            //pts.writeEndDivider(os);
        }


        // Write field
        OFstream(outputFile)() << tfield();
    }

    wroteGeom_ = true;
    return surfaceDir;
}


// Field writing methods
defineSurfaceWriterWriteFields(Foam::surfaceWriters::boundaryDataWriter);


// ************************************************************************* //
