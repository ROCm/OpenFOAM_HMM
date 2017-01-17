/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
     \\/     M anipulation  |
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
    foamVtkOutput

Description
    A collection of functions for writing vtk file content.

SourceFiles
    foamVtkOutput.C
    foamVtkOutputTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkOutput_H
#define foamVtkOutput_H

#include "floatScalar.H"
#include "volFields.H"
#include "foamVtkFormatter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class foamVtkOutput Declaration
\*---------------------------------------------------------------------------*/

class foamVtkOutput
{
    // Private Member Functions

        //- Disallow construction
        foamVtkOutput() = delete;

public:

    // Public typedefs

        //- Use UInt64 for header data
        typedef foamVtkFormatter::headerType headerType;


    // Forward declarations
    class legacy;


    // Static Members

        //- Write vtm datasets for specified files
        static Foam::label writeVtmFile
        (
            std::ostream& os,
            const UList<fileName>& files
        );


        //- Write a value component-wise.
        template<class Type>
        inline static void write(foamVtkFormatter&, const Type&);


        //- Write a list of values.
        //  The output does not include the payload size.
        template<class Type>
        static void writeList
        (
            foamVtkFormatter&,
            const UList<Type>&
        );


        //- Write a list of values via indirect addressing.
        //  The output does not include the payload size.
        template<class Type>
        static void writeList
        (
            foamVtkFormatter&,
            const UList<Type>&,
            const UList<label>& addressing
        );


        //- Write volField with cell values (including decomposed cells).
        //  The output includes the payload size and flush.
        template<class Type>
        static void writeField
        (
            foamVtkFormatter&,
            const GeometricField<Type, fvPatchField, volMesh>&,
            const UList<label>& superCells
        );

};


/*---------------------------------------------------------------------------*\
                    Class foamVtkOutput::legacy Declaration
\*---------------------------------------------------------------------------*/

//- Basic support for legacy files
class foamVtkOutput::legacy
{
    // Private Member Functions

        //- Disallow construction
        legacy() = delete;

public:

    // Static data members

        //- file extension for legacy files (vtk)
        static const Foam::word EXT;


    // Static Members

        //- Emit header for legacy file
        static std::ostream& writeHeader
        (
            std::ostream&,
            const std::string& title,
            const bool binary = false
        );


        //- Emit header for legacy CELL_DATA
        static std::ostream& writeCellDataHeader
        (
            std::ostream& os,
            const label nCells,
            const label nFields
        );


        //- Emit header for legacy POINT_DATA
        static std::ostream& writePointDataHeader
        (
            std::ostream& os,
            const label nPoints,
            const label nFields
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//- Template specialization for label
template<>
inline void Foam::foamVtkOutput::write<label>
(
    foamVtkFormatter& fmt,
    const label& val
)
{
    fmt.write(val);
}


//- Template specialization for float
template<>
inline void Foam::foamVtkOutput::write<float>
(
    foamVtkFormatter& fmt,
    const float& val
)
{
    fmt.write(val);
}


//- Template specialization for double
template<>
inline void Foam::foamVtkOutput::write<double>
(
    foamVtkFormatter& fmt,
    const double& val
)
{
    fmt.write(val);
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamVtkOutputTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
