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
    int32_t

Description
    32bit signed integer

SourceFiles
    int32.C
    int32IO.C

\*---------------------------------------------------------------------------*/

#ifndef int32_H
#define int32_H

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

//- Return a word representation of an int32
inline word name(const int32_t val)
{
    // No stripping required
    return word(std::to_string(val), false);
}


//- Return a word representation of an int32, using printf-style formatter.
//  The representation is not checked for valid word characters.
word name(const char* fmt, const int32_t val);


//- Return a word representation of an int32, using printf-style formatter.
//  The representation is not checked for valid word characters.
word name(const std::string& fmt, const int32_t val);


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Read int32_t from stream
int32_t readInt32(Istream& is);

//- Parse entire buffer as a int32_t, skipping leading/trailing whitespace.
//  \return Parsed value or FatalIOError on any problem
int32_t readInt32(const char* buf);

//- Parse entire string as a int32_t, skipping leading/trailing whitespace.
//  \return Parsed value or FatalIOError on any problem
inline int32_t readInt32(const std::string& str)
{
    return readInt32(str.c_str());
}

//- Read entire buffer as a int32_t, skipping leading/trailing whitespace.
//  \return True if successful.
bool readInt32(const char* buf, int32_t& val);

//- Read entire string as a int32_t, skipping leading/trailing whitespace.
//  \return True if successful.
inline bool readInt32(const std::string& str, int32_t& val)
{
    return readInt32(str.c_str(), val);
}

//- Same as readInt32
//  \return True if successful.
inline bool read(const char* buf, int32_t& val)
{
    return readInt32(buf, val);
}

//- Same as readInt32
//  \return True if successful.
inline bool read(const std::string& str, int32_t& val)
{
    return readInt32(str, val);
}


Istream& operator>>(Istream& is, int32_t& val);
Ostream& operator<<(Ostream& os, const int32_t val);

// On 32bit OSs long is not unambiguously int32_t (or int64_t) causing problems
// for IO operator resolution.
// This problem is avoided by explicitly defining the following operators:
#if WM_ARCH_OPTION == 32
    Istream& operator>>(Istream& is, long& val);
    Ostream& operator<<(Ostream& os, const long val);
#endif


//- Template specialization for pTraits<int32_t>
template<>
class pTraits<int32_t>
{
    int32_t p_;

public:

    //- Component type
    typedef int32_t cmptType;


    // Member constants

        //- Dimensionality of space
        static const direction dim = 3;

        //- Rank of int32_t is 0
        static const direction rank = 0;

        //- Number of components in int32_t is 1
        static const direction nComponents = 1;


    // Static data members

        static const char* const typeName;
        static const char* const componentNames[];
        static const int32_t zero;
        static const int32_t one;
        static const int32_t min;
        static const int32_t max;
        static const int32_t rootMax;
        static const int32_t rootMin;


    // Constructors

        //- Construct from primitive
        explicit pTraits(const int32_t& val);

        //- Construct from Istream
        pTraits(Istream& is);


    // Member Functions

        //- Access to the int32_t value
        operator int32_t() const
        {
            return p_;
        }

        //- Access to the int value
        operator int32_t&()
        {
            return p_;
        }
};


inline int32_t mag(const int32_t val)
{
    return ::abs(val);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
