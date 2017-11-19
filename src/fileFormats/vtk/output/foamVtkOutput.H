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
    Foam::vtk

Description
    Namespace for handling VTK output.
    Contains classes and functions for writing VTK file content.

Namespace
    Foam::vtk::legacy

Description
    Namespace for legacy VTK output constants and functions.

SourceFiles
    foamVtkOutput.C
    foamVtkOutputTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkOutput_H
#define foamVtkOutput_H

#include "autoPtr.H"
#include "Enum.H"
#include "foamVtkCore.H"
#include "foamVtkFormatter.H"
#include "floatScalar.H"
#include "IOstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  // Constants

  // General Functions

    //- Return a default asciiFormatter
    autoPtr<vtk::formatter> newFormatter(std::ostream& os);

    //- Return a new formatter based on the specified format type
    autoPtr<vtk::formatter> newFormatter
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
        vtk::formatter& fmt,
        const Type& val
    );


    //- Write a list of values.
    //  The output does not include the payload size.
    template<class Type>
    void writeList
    (
        vtk::formatter& fmt,
        const UList<Type>& lst
    );

    //- Write a list of values.
    //  The output does not include the payload size.
    template<class Type, unsigned Size>
    void writeList
    (
        vtk::formatter& fmt,
        const FixedList<Type, Size>& lst
    );


    //- Write a list of values via indirect addressing.
    //  The output does not include the payload size.
    template<class Type>
    void writeList
    (
        vtk::formatter& fmt,
        const UList<Type>& lst,
        const labelUList& addressing
    );


/*---------------------------------------------------------------------------*\
                               Namespace legacy
\*---------------------------------------------------------------------------*/

//- Some minimal additional support for writing legacy files
namespace legacy
{

  // Constants

    //- Strings corresponding to the (POLYDATA, UNSTRUCTURED_GRID) elements
    extern const Foam::Enum<vtk::fileTag> contentNames;

    //- Strings corresponding to the (CELL_DATA, POINT_DATA) elements
    extern const Foam::Enum<vtk::fileTag> dataTypeNames;


  // Functions

    //- Emit header for legacy file.
    //  Writes "ASCII" or "BINARY" depending on specified type.
    inline std::ostream& fileHeader
    (
        std::ostream& os,
        const std::string& title,
        const bool binary
    );

    //- Emit header for legacy file, with "ASCII" or "BINARY" depending on
    //  the formatter type.
    //  Includes "DATASET" with the specified dataset type.
    inline void fileHeader
    (
        vtk::formatter& fmt,
        const std::string& title,
        const vtk::fileTag& contentTypeTag
    );

    //- Emit header for legacy file, with "ASCII" or "BINARY" depending on
    //  the formatter type.
    //  If the contentType is non-empty, it is used for "DATASET" line.
    std::ostream& fileHeader
    (
        vtk::formatter& fmt,
        const std::string& title,
        const std::string& contentType
    );


    //- Emit header for POINTS (with trailing newline).
    inline void beginPoints(std::ostream& os, const label nPoints);

    //- Emit header for POLYGONS (with trailing newline).
    //  The nConnectivity is the sum of all connectivity points used,
    //  but \b without additional space for the size prefixes.
    //  The additional prefix sizes are added internally.
    inline void beginPolys
    (
        std::ostream& os,
        const label nPolys,
        const label nConnectivity
    );


    //- Use the enumerations vtk::fileTag::CELL_DATA, vtk::fileTag::POINT_DATA,
    //  to emit a legacy CELL_DATA, POINT_DATA element.
    //  The nEntries corresponds similarly to the number of cells or points,
    //  respectively.
    inline void dataHeader
    (
        std::ostream& os,
        const vtk::fileTag& dataTypeTag,
        const label nEntries,
        const label nFields
    );

    //- Start output of float field with the specified name.
    inline void floatField
    (
        std::ostream& os,
        const word& name,
        const int nCmpt,
        const label nEntries
    );

    //- Start output of int field with the specified name.
    inline void intField
    (
        std::ostream& os,
        const word& name,
        const int nCmpt,
        const label nEntries
    );

} // End namespace legacy


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "foamVtkOutputI.H"

#ifdef NoRepository
    #include "foamVtkOutputTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
