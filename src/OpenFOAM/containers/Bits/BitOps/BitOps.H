/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::BitOps

Description
    Various bit-wise operations, etc.

    The population count uses the Hamming weight
    (http://en.wikipedia.org/wiki/Hamming_weight).

\*---------------------------------------------------------------------------*/

#ifndef BitOps_H
#define BitOps_H

#include "label.H"
#include "Ostream.H"
#include <limits>
#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Namespace BitOps Declaration
\*---------------------------------------------------------------------------*/

namespace BitOps
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Count arbitrary number of bits (of an integral type)
template<class UIntType>
inline unsigned int bit_count(UIntType x)
{
    unsigned int n = 0u;

    for (; x; ++n) { x &= (x-1); }

    return n;
}


//- Count bits in a 32-bit value (Hamming weight method)
template<>
inline unsigned int bit_count(uint32_t x)
{
    x -= (x >> 1) & 0x55555555;
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333);

    return ((((x + (x >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24);
}


//- Count bits in a 64-bit value (Hamming weight method)
template<>
inline unsigned int bit_count(uint64_t x)
{
    x -= (x >> 1) & 0x5555555555555555;
    x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333);

    return unsigned
        ((((x + (x >> 4)) & 0x0F0F0F0F0F0F0F0F) * 0x0101010101010101) >> 56);
}


//- Repeat a value of the given BitWidth into the destination output type.
template<class UIntType, unsigned BitWidth>
inline UIntType repeat_value(unsigned val)
{
    static_assert
    (
        BitWidth && std::numeric_limits<UIntType>::digits >= BitWidth,
        "BitWidth too large for target output"
    );

    // How many fit into the target
    const unsigned nrepeat = (std::numeric_limits<UIntType>::digits / BitWidth);

    // Max value for a single element
    const unsigned mask = ((1u << BitWidth) - 1);

    // The first occurance
    UIntType fillval = ((val >= mask) ? mask : val);

    // Repeated
    for (unsigned i = 1; i < nrepeat; ++i)
    {
        fillval |= (fillval << BitWidth);
    }

    return fillval;
}


//- Print 0/1 bits in the (unsigned) integral type
template<class UIntType>
inline Ostream& print(Ostream& os, UIntType value, char off='0', char on='1')
{
    if (os.format() == IOstream::BINARY)
    {
        // Perhaps not the most sensible, but the only thing we currently have.
        os << label(value);
    }
    else
    {
        // Starting from most significant bit - makes for easy reading.
        for
        (
            unsigned test = (1u << (std::numeric_limits<UIntType>::digits-1));
            test;
            test >>= 1u
        )
        {
            os << ((value & test) ? on : off);
        }
    }

    return os;
}


//- An (unsigned) integral type adapter, for output of bit values
template<class UIntType>
struct bitInfo
{
    typedef UIntType value_type;

    value_type value;

    //- Null constructible as zero
    bitInfo() : value(0) {}

    //- Value construct
    explicit bitInfo(UIntType val) : value(val) {}

    //- Conversion to base type
    operator UIntType () const { return value; }

    //- Conversion to base type
    operator UIntType& () { return value; }
};


} // End namespace BitOps


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Print 0/1 bits of an (unsigned) integral type via an adapter
template<class UIntType>
inline Ostream& operator<<(Ostream& os, const BitOps::bitInfo<UIntType>& info)
{
    BitOps::print(os, info.value);
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
