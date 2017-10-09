/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

Primitive
    uint64_t

Description
    64bit unsigned integer

SourceFiles
    uint64.C
    uint64IO.C

\*---------------------------------------------------------------------------*/

#ifndef uint64_H
#define uint64_H

#include <cstdint>
#include <climits>
#include <cstdlib>

#include "word.H"
#include "pTraits.H"
#include "direction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Return a word representation of a uint64
inline word name(const uint64_t val)
{
    // No stripping required
    return word(std::to_string(val), false);
}


//- Return a word representation of a uint64_t, using printf-style formatter.
//  The representation is not checked for valid word characters.
word name(const char* fmt, const uint64_t);


//- Return a word representation of a uint64_t, using printf-style formatter.
//  The representation is not checked for valid word characters.
word name(const std::string& fmt, const uint64_t);


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Read uint64_t from stream.
uint64_t readUint64(Istream& is);

//- Parse entire buffer as uint64_t, skipping leading/trailing whitespace.
//  \return Parsed value or FatalIOError on any problem
uint64_t readUint64(const char* buf);

//- Parse entire string as uint64_t, skipping leading/trailing whitespace.
//  \return Parsed value or FatalIOError on any problem
inline uint64_t readUint64(const std::string& str)
{
    return readUint64(str.c_str());
}

//- Parse entire buffer as uint64_t, skipping leading/trailing whitespace.
//  \return True if successful.
bool readUint64(const char* buf, uint64_t& val);

//- Parse entire string as uint64_t, skipping leading/trailing whitespace.
//  \return True if successful.
inline bool readUint64(const std::string& str, uint64_t& val)
{
    return readUint64(str.c_str(), val);
}

//- Same as readUint64
//  \return True if successful.
inline bool read(const char* buf, uint64_t& val)
{
    return readUint64(buf, val);
}

//- Same as readUint64
//  \return True if successful.
inline bool read(const std::string& str, uint64_t& val)
{
    return readUint64(str, val);
}


Istream& operator>>(Istream& is, uint64_t& val);
Ostream& operator<<(Ostream& os, const uint64_t val);

//- Template specialization for pTraits<uint64_t>
template<>
class pTraits<uint64_t>
{
    uint64_t p_;

public:

    //- Component type
    typedef uint64_t cmptType;


    // Member constants

        //- Dimensionality of space
        static const direction dim = 3;

        //- Rank of uint64_t is 0
        static const direction rank = 0;

        //- Number of components in uint64_t is 1
        static const direction nComponents = 1;


    // Static data members

        static const char* const typeName;
        static const char* const componentNames[];
        static const uint64_t zero;
        static const uint64_t one;
        static const uint64_t min;
        static const uint64_t max;
        static const uint64_t rootMax;
        static const uint64_t rootMin;


    // Constructors

        //- Construct from primitive
        explicit pTraits(const uint64_t& val);

        //- Construct from Istream
        pTraits(Istream& is);


    // Member Functions

        //- Access to the uint64_t value
        operator uint64_t() const
        {
            return p_;
        }

        //- Access to the uint64_t value
        operator uint64_t&()
        {
            return p_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
