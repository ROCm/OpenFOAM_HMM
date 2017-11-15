/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::zero

Description
    A class representing the concept of 0 used to avoid unnecessary
    manipulations for objects that are known to be zero at compile-time.

SourceFiles
    zeroI.H

\*---------------------------------------------------------------------------*/

#ifndef zero_H
#define zero_H

#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class Istream;
class Ostream;

class zero;

/*---------------------------------------------------------------------------*\
                           Class zero Declaration
\*---------------------------------------------------------------------------*/

class zero
{
public:
    typedef zero value_type;

    // Forward declaration
    class null;


    // Constructors

        //- Null constructible
        zero()
        {}

        //- Construct null from Istream. Consumes no content.
        explicit zero(Istream&)
        {}


    // Member operators

        //- Return 0 for bool
        inline operator bool() const
        {
            return 0;
        }

        //- Return 0 for label
        inline operator label() const
        {
            return 0;
        }

        //- Return 0 for float
        inline operator float() const
        {
            return 0;
        }

        //- Return 0 for double
        inline operator double() const
        {
            return 0;
        }

};


/*---------------------------------------------------------------------------*\
                         Class zero::null Declaration
\*---------------------------------------------------------------------------*/

//- A zero class with a null output adapter.
class zero::null
:
    public zero
{
public:
    typedef null value_type;

    //- Null constructible
    null()
    {}

    //- Construct null from Istream without consuming any content.
    explicit null(Istream&)
    {}
};


// IOstream Operators

//- Read from Istream consumes no content.
inline Istream& operator>>(Istream& is, zero&)
{
    return is;
}

//- Write to Ostream emits no content.
inline Ostream& operator<<(Ostream& os, const zero::null&)
{
    return os;
}


// Global zero
static const zero Zero;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "zeroI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
