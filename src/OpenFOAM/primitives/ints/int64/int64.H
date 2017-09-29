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
    int64_t

Description
    64bit signed integer

SourceFiles
    int64.C
    int64IO.C

\*---------------------------------------------------------------------------*/

#ifndef int64_H
#define int64_H

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

//- Return a word representation of an int64
inline word name(const int64_t val)
{
    // No stripping required
    return word(std::to_string(val), false);
}


//- Return a word representation of an int64, using printf-style formatter.
//  The representation is not checked for valid word characters.
word name(const char* fmt, const int64_t);


//- Return a word representation of an int64, using printf-style formatter.
//  The representation is not checked for valid word characters.
word name(const std::string& fmt, const int64_t);


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Read int64_t from stream
int64_t readInt64(Istream& is);

//- Parse entire buffer as a int64_t, skipping leading/trailing whitespace.
//  \return Parsed value or FatalIOError on any problem
int64_t readInt64(const char* buf);

//- Parse entire string as a int64_t, skipping leading/trailing whitespace.
//  \return Parsed value or FatalIOError on any problem
inline int64_t readInt64(const std::string& str)
{
    return readInt64(str.c_str());
}

//- Read entire buffer as a int64_t, skipping leading/trailing whitespace.
//  \return True if successful.
bool readInt64(const char* buf, int64_t& val);

//- Read entire string as a int64_t, skipping leading/trailing whitespace.
//  \return True if successful.
inline bool readInt64(const std::string& str, int64_t& val)
{
    return readInt64(str.c_str(), val);
}

//- Same as readInt64
//  \return True if successful.
inline bool read(const char* buf, int64_t& val)
{
    return readInt64(buf, val);
}

//- Same as readInt64
//  \return True if successful.
inline bool read(const std::string& str, int64_t& val)
{
    return readInt64(str, val);
}


Istream& operator>>(Istream& is, int64_t& val);
Ostream& operator<<(Ostream& os, const int64_t val);

//- Template specialization for pTraits<int64_t>
template<>
class pTraits<int64_t>
{
    int64_t p_;

public:

    //- Component type
    typedef int64_t cmptType;


    // Member constants

        //- Dimensionality of space
        static const direction dim = 3;

        //- Rank of int64_t is 0
        static const direction rank = 0;

        //- Number of components in int64_t is 1
        static const direction nComponents = 1;


    // Static data members

        static const char* const typeName;
        static const char* const componentNames[];
        static const int64_t zero;
        static const int64_t one;
        static const int64_t min;
        static const int64_t max;
        static const int64_t rootMax;
        static const int64_t rootMin;


    // Constructors

        //- Construct from primitive
        explicit pTraits(const int64_t& val);

        //- Construct from Istream
        pTraits(Istream& is);


    // Member Functions

        //- Access to the int64_t value
        operator int64_t() const
        {
            return p_;
        }

        //- Access to the int value
        operator int64_t&()
        {
            return p_;
        }
};


inline int64_t mag(const int64_t val)
{
    return ::labs(val);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
