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

Class
    Foam::rawSurfaceWriter

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

/*---------------------------------------------------------------------------*\
                      Class rawSurfaceWriter Declaration
\*---------------------------------------------------------------------------*/

class rawSurfaceWriter
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
            const fileName& outputDir,      //!< output-dir
            const fileName& surfaceName,    //!< Name of surface
            const meshedSurf& surf,         //!< Surface geometry
            const word& fieldName,          //!< Name of field
            const Field<Type>& values,      //!< Field values to write
            const bool isNodeValues = false,//!< Values are per-vertex
            const bool verbose = false      //!< Additional verbosity
        ) const;


public:

    //- Runtime type information
    TypeName("raw");


    // Constructors

        //- Construct null
        rawSurfaceWriter();

        //- Construct with some output options
        rawSurfaceWriter(const dictionary& options);


    //- Destructor
    virtual ~rawSurfaceWriter() = default;


    // Member Functions

        //- Write single surface geometry to file.
        virtual fileName write
        (
            const fileName& outputDir,      //!< output-dir
            const fileName& surfaceName,    //!< Name of surface
            const meshedSurf& surf,         //!< Surface geometry
            const bool verbose = false      //!< Additional verbosity
        ) const; // override


        //- Write scalarField for a single surface to file.
        //  One value per face or vertex.
        virtual fileName write
        (
            const fileName& outputDir,      //!< output-dir
            const fileName& surfaceName,    //!< Name of surface
            const meshedSurf& surf,         //!< Surface geometry
            const word& fieldName,          //!< Name of field
            const Field<scalar>& values,    //!< Field values to write
            const bool isNodeValues = false,//!< Values are per-vertex
            const bool verbose = false      //!< Additional verbosity
        ) const; // override

        //- Write vectorField for a single surface to file.
        //  One value per face or vertex.
        virtual fileName write
        (
            const fileName& outputDir,      //!< output-dir
            const fileName& surfaceName,    //!< Name of surface
            const meshedSurf& surf,         //!< Surface geometry
            const word& fieldName,          //!< Name of field
            const Field<vector>& values,    //!< Field values to write
            const bool isNodeValues = false,//!< Values are per-vertex
            const bool verbose = false      //!< Additional verbosity
        ) const; // override

        //- Write sphericalTensorField for a single surface to file.
        //  One value per face or vertex.
        virtual fileName write
        (
            const fileName& outputDir,      //!< output-dir
            const fileName& surfaceName,    //!< Name of surface
            const meshedSurf& surf,         //!< Surface geometry
            const word& fieldName,          //!< Name of field
            const Field<sphericalTensor>& values, //!< Field values to write
            const bool isNodeValues = false,//!< Values are per-vertex
            const bool verbose = false      //!< Additional verbosity
        ) const; // override

        //- Write symmTensorField for a single surface to file.
        //  One value per face or vertex.
        virtual fileName write
        (
            const fileName& outputDir,      //!< output-dir
            const fileName& surfaceName,    //!< Name of surface
            const meshedSurf& surf,         //!< Surface geometry
            const word& fieldName,          //!< Name of field
            const Field<symmTensor>& values,//!< Field values to write
            const bool isNodeValues = false,//!< Values are per-vertex
            const bool verbose = false      //!< Additional verbosity
        ) const; // override

        //- Write tensorField for a single surface to file.
        //  One value per face or vertex.
        virtual fileName write
        (
            const fileName& outputDir,      //!< output-dir
            const fileName& surfaceName,    //!< Name of surface
            const meshedSurf& surf,         //!< Surface geometry
            const word& fieldName,          //!< Name of field
            const Field<tensor>& values,    //!< Field values to write
            const bool isNodeValues = false,//!< Values are per-vertex
            const bool verbose = false      //!< Additional verbosity
        ) const; // override
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
