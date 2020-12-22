/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
    Foam::surfaceWriters::nastranWriter

Description
    A surface writer for the Nastran file format - both surface mesh and fields

    The formatOptions for nastran:
    \table
        Property | Description                              | Required | Default
        fields   | field pairs for PLOAD2/PLOAD4            | yes   |
        format   | short / long / free                      | no    | long
        scale    | output geometry scaling                  | no    | 1
        fieldScale | output field scaling (dictionary)      | no    | empty
        commonGeometry | use separate geometry files        | no    | false
    \endtable

    For example,
    \verbatim
    formatOptions
    {
        nastran
        {
            // OpenFOAM field name to NASTRAN load types
            fields
            (
                (pMean  PLOAD2)
                (p      PLOAD4)
            );

            format  free;     // format type

            scale   1000;     // [m] -> [mm]
            fieldScale
            {
               "p.*"   0.01;  // [Pa] -> [mbar]
            }
        }
    }
    \endverbatim

    \section Output file locations

    The \c rootdir normally corresponds to something like
    \c postProcessing/\<name\>

    \subsection Geometry
    \verbatim
    rootdir
    `-- <time>
        |-- surfaceName0.{nas}
        `-- surfaceName1.{nas}
    \endverbatim

    \subsection Fields
    \verbatim
    rootdir
    `-- <time>
        `-- field0
        |   |-- surfaceName0.{bdf}
        |   `-- surfaceName1.{bdf}
        `-- field1
            |-- surfaceName0.{bdf}
            `-- surfaceName1.{bdf}
    \endverbatim

Note
    Output variable scaling does not apply to integer types such as Ids.
    Field pairs default to PLOAD2 for scalars and PLOAD4 for vectors etc.

SourceFiles
    nastranSurfaceWriter.C
    nastranSurfaceWriterImpl.C

\*---------------------------------------------------------------------------*/

#ifndef nastranSurfaceWriter_H
#define nastranSurfaceWriter_H

#include "surfaceWriter.H"
#include "NASCore.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{

/*---------------------------------------------------------------------------*\
                        Class nastranWriter Declaration
\*---------------------------------------------------------------------------*/

class nastranWriter
:
    public surfaceWriter
{
public:

    //- File field formats
    using fieldFormat = Foam::fileFormats::NASCore::fieldFormat;

    //- Output load format
    using loadFormat = Foam::fileFormats::NASCore::loadFormat;


private:

    // Private Data

        //- Field format (width and separator)
        fieldFormat writeFormat_;

        //- Mapping from field name to data format enumeration
        HashTable<loadFormat> fieldMap_;

        //- Use common geometry file
        bool commonGeometry_;

        //- Output geometry scaling
        const scalar geometryScale_;

        //- Output field scaling
        const dictionary fieldScale_;

        //- Separator used for free format
        word separator_;


    // Private Member Functions

        //- Write a coordinate
        void writeCoord
        (
            Ostream& os,
            const point& pt,
            const label pointI      //!< 0-based Point Id
        ) const;

        //- Write a face element (CTRIA3 or CQUAD4)
        void writeFace
        (
            Ostream& os,
            const word& faceType,
            const labelUList& facePts,
            const label elemId,     //!< 1-based Element Id
            const label propId      //!< 1-based Property Id
        ) const;

        //- Write the surface mesh geometry, tracking face decomposition
        //  Includes SHELL/MAT information
        //
        //  \param decompOffsets  begin/end offsets (size+1) into decompFaces
        //  \param decompFaces    Non tri/quad decomposed into triangles
        void writeGeometry
        (
            Ostream& os,
            const meshedSurf& surf,
            labelList& decompOffsets,
            DynamicList<face>& decompFaces
        ) const;

        //- Write the formatted keyword to the output stream
        Ostream& writeKeyword
        (
            Ostream& os,
            const word& keyword
        ) const;

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
            const label elemId      //!< 1-based Element Id
        ) const;


        //- Templated write operation
        template<class Type>
        fileName writeTemplate
        (
            const word& fieldName,          //!< Name of field
            const Field<Type>& localValues  //!< Local field values to write
        );


public:

    //- Declare type-name, virtual type (with debug switch)
    TypeNameNoDebug("nastran");


    // Constructors

        //- Default construct. Default SHORT format
        nastranWriter();

        //- Construct with some output options. Default LONG format
        explicit nastranWriter(const dictionary& options);

        //- Construct from components
        nastranWriter
        (
            const meshedSurf& surf,
            const fileName& outputPath,
            bool parallel = Pstream::parRun(),
            const dictionary& options = dictionary()
        );

        //- Construct from components
        nastranWriter
        (
            const pointField& points,
            const faceList& faces,
            const fileName& outputPath,
            bool parallel = Pstream::parRun(),
            const dictionary& options = dictionary()
        );


    //- Destructor
    virtual ~nastranWriter() = default;


    // Member Functions

        //- Format uses faceIds as part of its output
        virtual bool usesFaceIds() const // override
        {
            return true;
        }

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
