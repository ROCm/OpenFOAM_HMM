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

Class
    Foam::vtkPTraits

Description
    Names for VTK primitive types.

\*---------------------------------------------------------------------------*/

#ifndef foamVtkPTraits_H
#define foamVtkPTraits_H

#include <cstdint>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class endian;

/*---------------------------------------------------------------------------*\
                         Class vtkPTraits Declaration
\*---------------------------------------------------------------------------*/

template<class PrimitiveType>
class vtkPTraits
{
public:

    // Static data members

        static const char* const typeName;
};


template<>
const char* const vtkPTraits<uint8_t>::typeName;  // UInt8

template<>
const char* const vtkPTraits<int32_t>::typeName;  // Int32

template<>
const char* const vtkPTraits<uint32_t>::typeName; // UInt32

template<>
const char* const vtkPTraits<int64_t>::typeName;  // Int64

template<>
const char* const vtkPTraits<uint64_t>::typeName; // UInt64

template<>
const char* const vtkPTraits<float>::typeName;    // Float32

template<>
const char* const vtkPTraits<double>::typeName;   // Float64

template<>
const char* const vtkPTraits<Foam::endian>::typeName; // (Big|Little)Endian


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
