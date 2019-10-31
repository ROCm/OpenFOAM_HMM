/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2019 OpenCFD Ltd.
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

Description
    Help with architecture-specific aspects.

    Defines WM_BIG_ENDIAN or WM_LITTLE_ENDIAN as well as providing a
    few runtime methods. Primarily used as a namespace, but defined as
    a class to allow pTraits specialization.

SourceFiles
    endianI.H

\*---------------------------------------------------------------------------*/

#ifndef foamEndian_H  // Prefixed with 'foam' to avoid any name clashes
#define foamEndian_H

#include <cstdint>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef __BYTE_ORDER__
    // Clang, Gcc, Icc, Pgi
    #if (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
        #define WM_LITTLE_ENDIAN
    #elif (__BYTE_ORDER__ == __ORDER_BIG_ENDIAN__)
        #define WM_BIG_ENDIAN
    #else
        #error "__BYTE_ORDER__ is not BIG or LITTLE endian"
    #endif
#endif

// Could downgrade to a warning, but then user always needs runtime check.
#if !defined(WM_BIG_ENDIAN) && !defined(WM_LITTLE_ENDIAN)
    #error "Cannot determine BIG or LITTLE endian."
    #error "Please add to compilation options"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class endian Declaration
\*---------------------------------------------------------------------------*/

class endian
{
public:

    // Public Methods

        //- Runtime check for big endian.
        inline static bool isBig();

        //- Runtime check for little endian.
        inline static bool isLittle();

        //- Byte endian swapping for 32-bits
        inline static uint32_t swap32(uint32_t);

        //- Byte endian swapping for 64-bits
        inline static uint64_t swap64(uint64_t);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "endianI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
