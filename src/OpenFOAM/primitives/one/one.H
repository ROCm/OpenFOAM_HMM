/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    A class representing the concept of 1 (one) that can be used to avoid
    manipulating objects known to be \em one at compile-time.
    It is also used for tagged dispatch.

SourceFiles
    oneI.H

SeeAlso
    Foam::zero

\*---------------------------------------------------------------------------*/

#ifndef one_H
#define one_H

#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class one;
class Istream;
class Ostream;

/*---------------------------------------------------------------------------*\
                             Class one Declaration
\*---------------------------------------------------------------------------*/

class one
{
public:

    typedef one value_type;

    // Forward Declarations
    class minus;
    class null;

    //- Default construct
    constexpr one() noexcept {}

    //- Construct from Istream consumes no content.
    explicit constexpr one(Istream&) noexcept {}


    //- Return 1 for label
    constexpr operator label() const noexcept
    {
        return 1;
    }

    //- Return 1 for float
    constexpr operator float() const noexcept
    {
        return 1;
    }

    //- Return 1 for double
    constexpr operator double() const noexcept
    {
        return 1;
    }
};


/*---------------------------------------------------------------------------*\
                         Class one::minus Declaration
\*---------------------------------------------------------------------------*/

//- A class representing the concept of -1.
//  \note this class must never be derived from Foam::one, since this could
//  potentially mask its behaviour.
class one::minus
{
public:

    typedef minus value_type;

    //- Default construct
    constexpr minus() noexcept {}

    //- Construct from Istream consumes no content.
    explicit constexpr minus(Istream&) noexcept {}


    //- Return -1 for label
    constexpr operator label() const noexcept
    {
        return -1;
    }

    //- Return -1 for float
    constexpr operator float() const noexcept
    {
        return -1;
    }

    //- Return -1 for double
    constexpr operator double() const noexcept
    {
        return -1;
    }
};


/*---------------------------------------------------------------------------*\
                          Class one::null Declaration
\*---------------------------------------------------------------------------*/

//- A Foam::one class with a null output adapter.
class one::null
:
    public one
{
public:

    typedef null value_type;

    //- Default construct
    constexpr null() noexcept {}

    //- Construct from Istream consumes no content.
    explicit constexpr null(Istream&) noexcept {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Global one (1)
static constexpr const one One;

// IOstream Operators

//- Read from Istream consumes no content.
inline constexpr Istream& operator>>(Istream& is, one&) noexcept
{
    return is;
}

//- Read from Istream consumes no content.
inline constexpr Istream& operator>>(Istream& is, one::minus&) noexcept
{
    return is;
}

//- Write to Ostream emits no content.
inline constexpr Ostream& operator<<(Ostream& os, const one::null&) noexcept
{
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Global Operators, Functions
#include "oneI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
