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
    Hashing.C

\*---------------------------------------------------------------------------*/

#ifndef Hashing_H
#define Hashing_H

#include <string>
#include <climits>
#include <cstddef>
#include <stdint.h>    // C++0x uses <cstdint>

#include "label.H"
#include "uLabel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Macro defines

/**
 * @def bitRotateLeft(val, nBits)
 * Left-rotate a 32-bit value and carry by nBits
*/
#define bitRotateLeft(x, nBits)  (((x) << (nBits)) | ((x) >> (32 - (nBits))))

/**
 * @def bitMixer(a, b, c)
 * Mix three 32-bit values reversibly.
*/
#define bitMixer(a, b, c)                                                     \
    {                                                                         \
        a -= c; a ^= bitRotateLeft(c, 4); c += b;                             \
        b -= a; b ^= bitRotateLeft(a, 6); a += c;                             \
        c -= b; c ^= bitRotateLeft(b, 8); b += a;                             \
        a -= c; a ^= bitRotateLeft(c,16); c += b;                             \
        b -= a; b ^= bitRotateLeft(a,19); a += c;                             \
        c -= b; c ^= bitRotateLeft(b, 4); b += a;                             \
    }

/**
 * @def bitMixerFinal(a, b, c)
 * Final mixing of three 32-bit values (a,b,c) into c
*/
#define bitMixerFinal(a, b, c)                                                \
    {                                                                         \
        c ^= b; c -= bitRotateLeft(b, 14);                                    \
        a ^= c; a -= bitRotateLeft(c, 11);                                    \
        b ^= a; b -= bitRotateLeft(a, 25);                                    \
        c ^= b; c -= bitRotateLeft(b, 16);                                    \
        a ^= c; a -= bitRotateLeft(c, 4);                                     \
        b ^= a; b -= bitRotateLeft(a, 14);                                    \
        c ^= b; c -= bitRotateLeft(b, 24);                                    \
    }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace Hashing
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //- Bob Jenkins's 96-bit mixer hashing function (lookup3)
    //  @param[in] data - an array of uint32_t values
    //  @param[in] len  - the number of values (not bytes)
    //  @param[in] seed - the previous hash, or an arbitrary value
    uint32_t jenkins(const uint32_t*, size_t length, uint32_t seed = 0);

    //- Bob Jenkins's 96-bit mixer hashing function returning twin hash values
    //  @param[in] data - an array of uint32_t values
    //  @param[in] len  - the number of values (not bytes)
    //  @param[in] hash1 - the previous hash, or an arbitrary value
    //    on output, the primary hash value
    //  @param[in] hash1 - the previous hash, or an arbitrary value
    //    on output, the secondary hash value
    uint32_t jenkinsTwin
    (
        const uint32_t*,
        size_t len,
        uint32_t& hash1,
        uint32_t& hash2
    );

    //- Bob Jenkins's 96-bit mixer hashing function (lookup3)
    //  @param[in] data - a character stream
    //  @param[in] len  - the number of bytes
    //  @param[in] seed - the previous hash, or an arbitrary value
    uint32_t jenkins(const void* data, size_t len, uint32_t seed = 0);

    //- Hashing a C++ string with Jenkins' 96-bit mixer
    inline uint32_t jenkins(const std::string& str, uint32_t seed = 0)
    {
        return jenkins(str.data(), str.size(), seed);
    }

    //- Hash an unsigned long like two 32bit values
    inline uint32_t jenkins(const unsigned long& val)
    {
        uint32_t a = val >> 32;
        uint32_t b = val & 0xFFFFFFFF;
        uint32_t c = 0;

        bitMixerFinal(a,b,c);
        return c;
    }

    //- Hash a (64-bit) pointer like any long
    //  Treats bits as zero if the pointer is actually 32bit
    inline uint32_t pointerHash(const void* const& p)
    {
        return jenkins(long(p));
    }

#if UINT_MAX > WANTEDURANGE

    //- Hash a uLabel
    inline unsigned intHash(const uLabel& val1)
    {
        uint32_t a = val1;
        uint32_t b = 0;
        uint32_t c = 0;

        bitMixerFinal(a,b,c);
        return c;
    }

    //- Hash a uLabel pair
    inline unsigned intHash
    (
        const uLabel& val1,
        const uLabel& val2
    )
    {
        uint32_t a = val1;
        uint32_t b = val2;
        uint32_t c = 0;

        bitMixerFinal(a,b,c);
        return c;
    }

    //- Hash a uLabel triplet
    inline unsigned intHash
    (
        const uLabel& val1,
        const uLabel& val2,
        const uLabel& val3
    )
    {
        uint32_t a = val1;
        uint32_t b = val2;
        uint32_t c = val3;

        bitMixer(a,b,c);

        bitMixerFinal(a,b,c);
        return c;
    }

#else

    //- Hash a uLabel
    inline unsigned intHash
    (
        const uLabel& val1
    )
    {
        uint32_t a = val1 >> 32;
        uint32_t b = val1 & 0xFFFFFFFF;
        uint32_t c = 0;

        bitMixerFinal(a,b,c);
        return c;
    }

    //- Hash a uLabel pair
    inline unsigned intHash
    (
        const uLabel& val1,
        const uLabel& val2
    )
    {
        uint32_t a = val1 >> 32;
        uint32_t b = val1 & 0xFFFFFFFF;
        uint32_t c = val2 >> 32;

        bitMixer(a,b,c);

        a += val2 & 0xFFFFFFFF;

        bitMixerFinal(a,b,c);
        return c;
    }

    //- Hash a uLabel triplet
    inline unsigned intHash
    (
        const uLabel& val1,
        const uLabel& val2,
        const uLabel& val3
    )
    {
        uint32_t a = val1 >> 32;
        uint32_t b = val1 & 0xFFFFFFFF;
        uint32_t c = val2 >> 32;

        bitMixer(a,b,c);

        a += val2 & 0xFFFFFFFF;
        b += val3 >> 32;
        c += val3 & 0xFFFFFFFF;

        bitMixerFinal(a,b,c);
        return c;
    }

#endif

    //- Hash a label as uLabel
    inline unsigned intHash(const label& val)
    {
        return intHash(static_cast<uLabel>(val));
    }

    //- Hash a label pair
    inline unsigned intHash
    (
        const label& val1,
        const label& val2
    )
    {
        return intHash
        (
            static_cast<uLabel>(val1),
            static_cast<uLabel>(val2)
        );
   }

    //- Hash a label triplet
    inline unsigned intHash
    (
        const label& val1,
        const label& val2,
        const label& val3
    )
    {
        return intHash
        (
            static_cast<uLabel>(val1),
            static_cast<uLabel>(val2),
            static_cast<uLabel>(val3)
        );
   }

#if 0
    inline unsigned intHashCommutative(const label s1, const label s2)
    {
        if (s1 < s2)
        {
            intHash(s1, s2);
        }
        else
        {
            intHash(s2, s1);
        }
    }

#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Hashing

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
