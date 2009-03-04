/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Namespace
    Foam::Hashing

Description
    Misc. hashing functions, mostly from Bob Jenkins.

    The Jenkins hashing function(s) is similar in speed to Paul Hsieh's
    SuperFast hash, but is public domain, supports incremental hashing
    and has been reported to have better characteristics.
    It is also what postgresql seems to be using.

SeeAlso
    http://burtleburtle.net/bob/c/lookup3.c

SourceFiles
    Hasher.C

\*---------------------------------------------------------------------------*/

#ifndef Hasher_H
#define Hasher_H

#include <climits>
#include <cstddef>
#include <stdint.h>    // C++0x uses <cstdint>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //- Bob Jenkins's 96-bit mixer hashing function (lookup3)
    //  @param[in] data - a character stream
    //  @param[in] len  - the number of bytes
    //  @param[in] seed - the previous hash, or an arbitrary value
    unsigned Hasher(const void* data, size_t len, unsigned seed = 0);

    //- An optimized version of Hasher
    //  @param[in] data - an array of uint32_t values
    //  @param[in] len  - the number of values (not bytes)
    //  @param[in] seed - the previous hash, or an arbitrary value
    unsigned HasherSingle(const uint32_t*, size_t length, unsigned seed = 0);

    //- An optimized version of Hasher, returning dual hash values
    //  @param[in] data - an array of uint32_t values
    //  @param[in] len  - the number of values (not bytes)
    //  @param[in] hash1 - the previous hash, or an arbitrary value
    //    on output, the primary hash value
    //  @param[in] hash1 - the previous hash, or an arbitrary value
    //    on output, the secondary hash value
    unsigned HasherDual
    (
        const uint32_t*,
        size_t len,
        unsigned& hash1,
        unsigned& hash2
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
