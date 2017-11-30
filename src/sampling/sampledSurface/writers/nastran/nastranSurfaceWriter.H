/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2017 OpenCFD Ltd.
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
    Foam::nastranSurfaceWriter

Description
    A surface writer for the Nastran file format - both surface mesh and fields

    \verbatim
    formatOptions
    {
        nastran
        {
            // From OpenFOAM field name to NASTRAN field name
            fields
            (
                (pMean  PLOAD2)
                (p      PLOAD4)
            );
            // Optional scale
            scale   2.0;
            // Optional format
            format  free;    // short, long, free
        }
    };
    \endverbatim

SourceFiles
    nastranSurfaceWriter.C
    nastranSurfaceWriterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef nastranSurfaceWriter_H
#define nastranSurfaceWriter_H

#include "surfaceWriter.H"
#include "NASCore.H"
#include "OFstream.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class nastranSurfaceWriter Declaration
\*---------------------------------------------------------------------------*/

class nastranSurfaceWriter
:
    public surfaceWriter
{
public:

    //- File field formats
    using fieldFormat = Foam::fileFormats::NASCore::fieldFormat;

    //- Output load format
    enum loadFormat
    {
        PLOAD2,
        PLOAD4
    };


private:

    // Private data

        static const Enum<loadFormat> loadFormatNames_;

        //- Field format (width and separator)
        fieldFormat writeFormat_;

        //- Mapping from field name to data format enumeration
        HashTable<loadFormat> fieldMap_;

        //- Scale to apply to values (default = 1.0)
        scalar scale_;

        //- Separator used for free format
        word separator_;


    // Private Member Functions

        //- Write a coordinate
        void writeCoord
        (
            Ostream& os,
            const point& pt,
            const label pointI  //!< 0-based Point Id
        ) const;

        //- Write a face element (CTRIA3 or CQUAD4)
        void writeFace
        (
            Ostream& os,
            const word& faceType,
            const labelUList& facePts,
            const label EID,    //!< 1-based Element Id
            const label PID     //!< 1-based Property Id
        ) const;

        //- Main driver to write the surface mesh geometry
        void writeGeometry
        (
            Ostream& os,
            const meshedSurf& surf,
            List<DynamicList<face>>& decomposedFaces
        ) const;

        //- Write the formatted keyword to the output stream
        Ostream& writeKeyword
        (
            Ostream& os,
            const word& keyword
        ) const;

        //- Write the footer information
        Ostream& writeFooter(Ostream& os, const meshedSurf& surf) const;

        //- Write a formatted value to the output stream
        template<class Type>
        Ostream& writeValue(Ostream& os, const Type& value) const;

        //- Write a face-based value
        template<class Type>
        Ostream& writeFaceValue
        (
            Ostream& os,
            const loadFormat format,
            const Type& value,
            const label EID     //!< 1-based Element Id
        ) const;


        //- Templated write operation
        template<class Type>
        fileName writeTemplate
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const meshedSurf& surf,
            const word& fieldName,
            const Field<Type>& values,
            const bool isNodeValues,
            const bool verbose
        ) const;


public:

    //- Runtime type information
    TypeName("nastran");


    // Constructors

        //- Construct null
        nastranSurfaceWriter();

        //- Construct with some output options
        nastranSurfaceWriter(const dictionary& options);


    //- Destructor
    virtual ~nastranSurfaceWriter() = default;


    // Member Functions

        //- True if the surface format supports geometry in a separate file.
        //  False if geometry and field must be in a single file
        virtual bool separateGeometry() const
        {
            return false;
        }

        //- Write single surface geometry to file.
        virtual fileName write
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const meshedSurf& surf,
            const bool verbose = false
        ) const; // override

        //- Write scalarField for a single surface to file.
        //  One value per face or vertex (isNodeValues = true)
        virtual fileName write
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const meshedSurf& surf,
            const word& fieldName,
            const Field<scalar>& values,
            const bool isNodeValues,
            const bool verbose = false
        ) const; // override

        //- Write vectorField for a single surface to file.
        //  One value per face or vertex (isNodeValues = true)
        virtual fileName write
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const meshedSurf& surf,
            const word& fieldName,
            const Field<vector>& values,
            const bool isNodeValues,
            const bool verbose = false
        ) const; // override

        //- Write sphericalTensorField for a single surface to file.
        //  One value per face or vertex (isNodeValues = true)
        virtual fileName write
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const meshedSurf& surf,
            const word& fieldName,
            const Field<sphericalTensor>& values,
            const bool isNodeValues,
            const bool verbose = false
        ) const; // override

        //- Write symmTensorField for a single surface to file.
        //  One value per face or vertex (isNodeValues = true)
        virtual fileName write
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const meshedSurf& surf,
            const word& fieldName,
            const Field<symmTensor>& values,
            const bool isNodeValues,
            const bool verbose = false
        ) const; // override

        //- Write tensorField for a single surface to file.
        //  One value per face or vertex (isNodeValues = true)
        virtual fileName write
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const meshedSurf& surf,
            const word& fieldName,
            const Field<tensor>& values,
            const bool isNodeValues,
            const bool verbose = false
        ) const; // override
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nastranSurfaceWriterTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
