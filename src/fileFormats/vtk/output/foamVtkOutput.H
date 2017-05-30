/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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

Namespace
    Foam::foamVtkOutput

Description
    Namespace for handling VTK output.
    Contains classes and functions for writing VTK file content.

Namespace
    Foam::foamVtkOutput::legacy

Description
    Namespace for legacy VTK output functions.

SourceFiles
    foamVtkOutput.C
    foamVtkOutputTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkOutput_H
#define foamVtkOutput_H

#include "autoPtr.H"
#include "Enum.H"
#include "foamVtkOutputTypes.H"
#include "foamVtkFormatter.H"
#include "floatScalar.H"
#include "IOstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace foamVtkOutput
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  // Typedefs

    //- Use UInt64 for header data
    typedef formatter::headerType headerType;


  // General Functions

    //- Return a default asciiFormatter
    autoPtr<foamVtkOutput::formatter> newFormatter(std::ostream& os);


    //- Return a new formatter based on the specified format type
    autoPtr<foamVtkOutput::formatter> newFormatter
    (
        std::ostream& os,
        const enum formatType fmtType,
        unsigned prec = IOstream::defaultPrecision()
    );


    //- Write vtm datasets for specified files
    label writeVtmFile(std::ostream& os, const UList<fileName>& files);


    //- Write a value component-wise.
    template<class Type>
    inline void write
    (
        foamVtkOutput::formatter& fmt,
        const Type& val
    );


    //- Write a list of values.
    //  The output does not include the payload size.
    template<class Type>
    void writeList
    (
        foamVtkOutput::formatter& fmt,
        const UList<Type>& lst
    );


    //- Write a list of values via indirect addressing.
    //  The output does not include the payload size.
    template<class Type>
    void writeList
    (
        foamVtkOutput::formatter& fmt,
        const UList<Type>& lst,
        const UList<label>& addressing
    );


/*---------------------------------------------------------------------------*\
                               Namespace legacy
\*---------------------------------------------------------------------------*/

//- Some minimal additional support for writing legacy files
namespace legacy
{
    //- Some standard text elements for legacy vtk files
    enum class textElemType
    {
        POINTS,             //!< "POINTS"
        CELLS,              //!< "CELLS"
        POLYS,              //!< "POLYGONS"
        VERTS,              //!< "VERTICES"
        POLY_DATA,          //!< "POLYDATA"
        UNSTRUCTURED_GRID,  //!< "UNSTRUCTURED_GRID"
    };

    //- Strings corresponding to the elements
    extern const Foam::Enum<textElemType> textNames;

    //- Strings corresponding to the elements
    extern const Foam::Enum<vtkFileTag>   dataTypeNames;


  // Functions

    //- Emit header for legacy file.
    //  Writes "ASCII" or "BINARY" depending on specified type.
    std::ostream& fileHeader
    (
        std::ostream& os,
        const std::string& title,
        const bool binary = false
    );

    //- Emit header for legacy file.
    //  Writes "ASCII" or "BINARY" depending on the formatter type.
    std::ostream& fileHeader
    (
        foamVtkOutput::formatter& fmt,
        const std::string& title
    );


    //- Emit header for legacy CELL_DATA or POINT_DATA, corresponding to the
    //  enumeration textElemType::CELLS or textElemType::POINTS, respectively.
    //  The nEntries corresponds similarly to the number of cells or points,
    //  respectively.
    std::ostream& dataHeader
    (
        std::ostream& os,
        const vtkFileTag& dataTypeTag,
        const label nEntries,
        const label nFields
    );

    //- Emit header for legacy CELL_DATA.
    //  The nEntries should normally correspond to the number of cells.
    std::ostream& cellDataHeader
    (
        std::ostream& os,
        const label nEntries,
        const label nFields
    );

    //- Emit header for legacy POINT_DATA
    //  The nEntries should normally correspond to the number of points.
    std::ostream& pointDataHeader
    (
        std::ostream& os,
        const label nEntries,
        const label nFields
    );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Template specialization for label
template<>
inline void write<label>(foamVtkOutput::formatter& fmt, const label& val)
{
    fmt.write(val);
}


//- Template specialization for float
template<>
inline void write<float>(foamVtkOutput::formatter& fmt, const float& val)
{
    fmt.write(val);
}


//- Template specialization for double
template<>
inline void write<double>(foamVtkOutput::formatter& fmt, const double& val)
{
    fmt.write(val);
}

} // End namespace foamVtkOutput


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamVtkOutputTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
