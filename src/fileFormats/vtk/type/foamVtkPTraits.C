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

\*---------------------------------------------------------------------------*/

#include "foamVtkPTraits.H"
#include "endian.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const
Foam::foamVtkPTraits<uint8_t>::typeName = "UInt8";

template<>
const char * const
Foam::foamVtkPTraits<int32_t>::typeName = "Int32";

template<>
const char * const
Foam::foamVtkPTraits<uint32_t>::typeName = "UInt32";

template<>
const char * const
Foam::foamVtkPTraits<int64_t>::typeName = "Int64";

template<>
const char * const
Foam::foamVtkPTraits<uint64_t>::typeName = "UInt64";

template<>
const char * const
Foam::foamVtkPTraits<float>::typeName = "Float32";

template<>
const char * const
Foam::foamVtkPTraits<double>::typeName = "Float64";

#ifdef WM_LITTLE_ENDIAN
template<>
const char* const
Foam::foamVtkPTraits<Foam::endian>::typeName = "LittleEndian";
#else
template<>
const char* const
Foam::foamVtkPTraits<Foam::endian>::typeName = "BigEndian";
#endif


// ************************************************************************* //
