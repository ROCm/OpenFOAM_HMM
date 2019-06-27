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

Class
    Foam::surfaceWriters::rawWriter

Description
    A surfaceWriter for raw output.

    The formatOptions for raw:
    \table
        Property     | Description                          | Required | Default
        compression  | on / off                             | no    | off
    \endtable

    For example,
    \verbatim
    formatOptions
    {
        raw
        {
            compression on;
        }
    }
    \endverbatim

    \heading Output file locations

    The \c rootdir normally corresponds to something like
    \c postProcessing/\<name\>

    \subheading Geometry
    \verbatim
    rootdir
    `-- timeName
        `-- surfaceName.{raw}
    \endverbatim

    \subheading Fields
    \verbatim
    rootdir
    `-- timeName
        |-- <field0>_surfaceName.{raw}
        `-- <field1>_surfaceName.{raw}
    \endverbatim

SourceFiles
    rawSurfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef rawSurfaceWriter_H
#define rawSurfaceWriter_H

#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{

/*---------------------------------------------------------------------------*\
                          Class rawWriter Declaration
\*---------------------------------------------------------------------------*/

class rawWriter
:
    public surfaceWriter
{
    // Private data

        //- Output compression (default: uncompressed)
        IOstream::compressionType writeCompression_;


    // Private Member Functions

        //- Templated write operation
        template<class Type>
        fileName writeTemplate
        (
            const word& fieldName,          //!< Name of field
            const Field<Type>& localValues  //!< Local field values to write
        );


public:

    //- Runtime type information
    TypeNameNoDebug("raw");


    // Constructors

        //- Construct null
        rawWriter();

        //- Construct with some output options
        explicit rawWriter(const dictionary& options);

        //- Construct from components
        rawWriter
        (
            const meshedSurf& surf,
            const fileName& outputPath,
            bool parallel = Pstream::parRun(),
            const dictionary& options = dictionary()
        );

        //- Construct from components
        rawWriter
        (
            const pointField& points,
            const faceList& faces,
            const fileName& outputPath,
            bool parallel = Pstream::parRun(),
            const dictionary& options = dictionary()
        );


    //- Destructor
    virtual ~rawWriter() = default;


    // Member Functions

        //- Write surface geometry to file.
        virtual fileName write(); // override

        declareSurfaceWriterWriteMethod(label);
        declareSurfaceWriterWriteMethod(scalar);
        declareSurfaceWriterWriteMethod(vector);
        declareSurfaceWriterWriteMethod(sphericalTensor);
        declareSurfaceWriterWriteMethod(symmTensor);
        declareSurfaceWriterWriteMethod(tensor);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceWriters
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
