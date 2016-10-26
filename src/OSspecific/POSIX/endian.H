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

Description
    Help with architecture-specific aspects.

    Defines WM_BIG_ENDIAN or WM_LITTLE_ENDIAN as well as providing a
    few runtime methods. Primarily used as a namespace, but provided
    as a class for possible future expansion.

SourceFiles
    endianI.H

\*---------------------------------------------------------------------------*/

#ifndef foamEndian_H  // prefixed with 'foam' to avoid potential collisions
#define foamEndian_H

#include <cstdint>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef __GNUC__
    #include <endian.h>
#elif defined __mips
    #include <standards.h>
    #include <sys/endian.h>
#endif

#ifdef __BYTE_ORDER
    #if (__BYTE_ORDER == __LITTLE_ENDIAN)
        #define WM_LITTLE_ENDIAN
    #elif (__BYTE_ORDER == __BIG_ENDIAN)
        #define WM_BIG_ENDIAN
    #else
        #error "__BYTE_ORDER defined, but not __LITTLE_ENDIAN or __BIG_ENDIAN."
    #endif
#elif defined(__LITTLE_ENDIAN)  \
   || defined(_LITTLE_ENDIAN)   \
   || defined(LITTLE_ENDIAN)
   #define WM_LITTLE_ENDIAN
#elif defined(__BIG_ENDIAN)     \
   || defined(_BIG_ENDIAN)      \
   || defined(BIG_ENDIAN)
   #define WM_BIG_ENDIAN
#endif

// Special handling for OS-X
#ifdef __DARWIN_BYTE_ORDER
    #if (__DARWIN_BYTE_ORDER == __DARWIN_BIG_ENDIAN)
        #define WM_BIG_ENDIAN
        #undef  WM_LITTLE_ENDIAN
    #elif (__DARWIN_BYTE_ORDER == __DARWIN_LITTLE_ENDIAN)
        #define WM_LITTLE_ENDIAN
        #undef  WM_BIG_ENDIAN
    #endif
#endif

// Could also downgrade to a warning, but then user always needs runtime check.
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
    // Private Member Functions

        //- Disallow default bitwise copy construct
        endian(const endian&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const endian&) = delete;

public:

    // Public data

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
