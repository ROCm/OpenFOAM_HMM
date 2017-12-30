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
    Foam::one

Description
    A class representing the concept of 1 (or 1.0) used to avoid unnecessary
    manipulations for objects that are known to be one at compile-time.

SourceFiles
    oneI.H

\*---------------------------------------------------------------------------*/

#ifndef one_H
#define one_H

#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes, friend functions and operators
class Istream;
class Ostream;
class one;

/*---------------------------------------------------------------------------*\
                           Class one Declaration
\*---------------------------------------------------------------------------*/

class one
{
public:
    typedef one value_type;

    // Forward declaration
    class null;


    // Constructors

        //- Null constructible
        one()
        {}

        //- Construct null from Istream. Consumes no content.
        explicit one(Istream&)
        {}


    // Member operators

        //- Return 1 for label
        inline operator label() const
        {
            return 1;
        }

        //- Return 1 for float
        inline operator float() const
        {
            return 1;
        }

        //- Return 1 for double
        inline operator double() const
        {
            return 1;
        }
};


/*---------------------------------------------------------------------------*\
                          Class one::null Declaration
\*---------------------------------------------------------------------------*/

//- A one class with a null output adapter.
class one::null
:
    public one
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
inline Istream& operator>>(Istream& is, one&)
{
    return is;
}

//- Write to Ostream emits no content.
inline Ostream& operator<<(Ostream& os, const one::null&)
{
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "oneI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
