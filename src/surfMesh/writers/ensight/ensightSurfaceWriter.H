/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2011, 2015-2019 OpenCFD Ltd.
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

Class
    Foam::surfaceWriters::ensightWriter

Description
    A surfaceWriter for Ensight format.

    \verbatim
    formatOptions
    {
        ensight
        {
            format          ascii;
            collateTimes    true;
        }
    }
    \endverbatim

    Format options:
    \table
        Property | Description                              | Required | Default
        format   | ascii/binary                             | no  | ascii
        collateTimes | use common geometry for times        | no  | true
    \endtable

SourceFiles
    ensightSurfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef ensightSurfaceWriter_H
#define ensightSurfaceWriter_H

#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{

/*---------------------------------------------------------------------------*\
                        Class ensightWriter Declaration
\*---------------------------------------------------------------------------*/

class ensightWriter
:
    public surfaceWriter
{
    // Private Data

        //- Write option (default: IOstream::ASCII)
        IOstream::streamFormat writeFormat_;

        //- Collate times (default: true)
        bool collateTimes_;


    // Private Member Functions

        //- Print time-set for ensight case file
        static void printTimeset
        (
            OSstream& os,
            const label ts,
            const scalar& timeValue
        );

        //- Print time-set for ensight case file
        static void printTimeset
        (
            OSstream& os,
            const label ts,
            const UList<scalar>& times
        );


        //- Write geometry
        fileName writeCollated();

        //- Write geometry
        fileName writeUncollated();

        //- Templated write operation - one file per timestep
        template<class Type>
        fileName writeCollated
        (
            const word& fieldName,          //!< Name of field
            const Field<Type>& localValues  //!< Local field values to write
        );

        //- Templated write operation - all time steps in single file
        template<class Type>
        fileName writeUncollated
        (
            const word& fieldName,          //!< Name of field
            const Field<Type>& localValues  //!< Local field values to write
        );

        //- Templated write operation
        template<class Type>
        fileName writeTemplate
        (
            const word& fieldName,          //!< Name of field
            const Field<Type>& localValues  //!< Local field values to write
        );

public:

    //- Runtime type information
    TypeNameNoDebug("ensight");


    // Constructors

        //- Construct null
        ensightWriter();

        //- Construct with some output options
        explicit ensightWriter(const dictionary& options);

        //- Construct from components
        ensightWriter
        (
            const meshedSurf& surf,
            const fileName& outputPath,
            bool parallel = Pstream::parRun(),
            const dictionary& options = dictionary()
        );

        //- Construct from components
        ensightWriter
        (
            const pointField& points,
            const faceList& faces,
            const fileName& outputPath,
            bool parallel = Pstream::parRun(),
            const dictionary& options = dictionary()
        );


    //- Destructor
    virtual ~ensightWriter() = default;


    // Member Functions

        //- True if the surface format supports geometry in a separate file.
        //  False if geometry and field must be in a single file
        virtual bool separateGeometry() const;

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
