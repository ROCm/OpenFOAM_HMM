/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    Foam::Hash

Description
    Hash function class.
    The default definition is for primitives.
    Non-primitives used to hash entries on hash tables will likely need
    a specialized version.

Note
    The second template parameter (bool) is used for SFINAE overloading,

\*---------------------------------------------------------------------------*/

#ifndef Hash_H
#define Hash_H

#include "Hasher.H"
#include <cstdint>
#include <string>
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class Hash Declaration
\*---------------------------------------------------------------------------*/

template<class T, class SFINAEType=bool>
struct Hash
{
    unsigned operator()(const T& obj, unsigned seed=0) const
    {
        return Foam::Hasher(&obj, sizeof(T), seed);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Specialization for trivial (integer) types

#undef  FOAM_INTHASHER
#define FOAM_INTHASHER(IntType)                                               \
    /*! \brief Hashing specialization for IntType */                          \
    /*! Unseeded (single value) uses natural order, otherwise incremental */  \
    template<> struct Hash<IntType>                                           \
    {                                                                         \
        unsigned operator()(const IntType val) const                          \
        {                                                                     \
            return static_cast<unsigned>(val);                                \
        }                                                                     \
        unsigned operator()(const IntType val, unsigned seed) const           \
        {                                                                     \
            return Foam::Hasher(&val, sizeof(IntType), seed);                 \
        }                                                                     \
    }

FOAM_INTHASHER(bool);
FOAM_INTHASHER(char);
FOAM_INTHASHER(int32_t);
FOAM_INTHASHER(int64_t);
FOAM_INTHASHER(uint32_t);

#undef FOAM_INTHASHER

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Hashing specialization for nullptr. Always 0
template<>
struct Hash<std::nullptr_t>
{
    unsigned operator()(std::nullptr_t, unsigned seed=0) const noexcept
    {
        return seed;
    }
};

//- Hashing specialization for pointers, interpret pointer as a integer type
template<>
struct Hash<void*>
{
    unsigned operator()(const void* const ptr, unsigned seed=0) const
    {
        const uintptr_t addr = uintptr_t(ptr);
        return Foam::Hasher(&addr, sizeof(addr), seed);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Hashing partial specialization for derived string types
template<class StringType>
struct Hash
<
    StringType,
    typename std::enable_if
    <std::is_base_of<std::string, StringType>::value, bool>::type
>
{
    unsigned operator()(const std::string& str, unsigned seed=0) const
    {
        return Foam::Hasher(str.data(), str.length(), seed);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
