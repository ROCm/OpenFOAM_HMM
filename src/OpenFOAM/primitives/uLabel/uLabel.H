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

Typedef
    Foam::uLabel

Description
    A uLabel is an unsigned label. See label.H.

\*---------------------------------------------------------------------------*/

#ifndef uLabel_H
#define uLabel_H

#include <climits>
#include <cstdlib>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#if FOAM_LABEL64
#    define WANTEDURANGE 18000000000000000000
#else
#    define WANTEDURANGE 4000000000
#endif


#if UINT_MAX > WANTEDURANGE

// Define uLabel as an unsigned int

#include "uint.H"

namespace Foam
{
    typedef unsigned int uLabel;

    static const uLabel uLabelMax = UINT_MAX;

    inline uLabel readULabel(Istream& is)
    {
        return readUint(is);
    }

} // End namespace Foam


#elif ULONG_MAX > WANTEDURANGE

// Define uLabel as an unsigned long

#include "uint.H"
#include "ulong.H"

namespace Foam
{
    typedef unsigned long uLabel;

    static const uLabel uLabelMax = LONG_MAX;

    inline uLabel readULabel(Istream& is)
    {
        return readUlong(is);
    }

} // End namespace Foam


#elif ULLONG_MAX > WANTEDURANGE

// Define uLabel as an unsigned long long

#error "Not implemented yet"

#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pTraits.H"
#include "direction.H"

namespace Foam
{

// template specialisation for pTraits<label>
template<>
class pTraits<uLabel>
{
    uLabel p_;

public:

    //- Component type
    typedef uLabel cmptType;

    // Member constants

        enum
        {
            dim = 3,         // Dimensionality of space
            rank = 0,        // Rank of uLabel is 0
            nComponents = 1  // Number of components in uLabel is 1
        };

    // Static data members

        static const char* const typeName;
        static const char* componentNames[];
        static const uLabel zero;
        static const uLabel one;

    // Constructors

        //- Construct from Istream
        pTraits(Istream&);


    // Member Functions

        operator uLabel() const
        {
            return p_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline uLabel& setComponent(uLabel& l, const direction)
{
    return l;
}

inline uLabel component(const uLabel l, const direction)
{
    return l;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
