/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
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
    A class representing the concept of 0 (zero), which can be used to avoid
    manipulating objects that are known to be \em zero at compile-time.

SourceFiles
    zero.C
    zeroI.H

SeeAlso
    Foam::one

\*---------------------------------------------------------------------------*/

#ifndef zero_H
#define zero_H

#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class zero;
class Istream;
class Ostream;

/*---------------------------------------------------------------------------*\
                            Class zero Declaration
\*---------------------------------------------------------------------------*/

class zero
{
public:
    typedef zero value_type;

    // Forward declarations
    class null;

    //- Null constructible
    constexpr zero() noexcept {}

    //- Construct from Istream consumes no content.
    explicit constexpr zero(Istream&) noexcept {}


    //- Return false (0) for bool
    inline constexpr operator bool() const noexcept
    {
        return false;
    }

    //- Return 0 for label
    inline constexpr operator label() const noexcept
    {
        return 0;
    }

    //- Return 0 for float
    inline constexpr operator float() const noexcept
    {
        return 0;
    }

    //- Return 0 for double
    inline constexpr operator double() const noexcept
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

    //- A static zero::null for dereferencing as a dummy element
    static null dummy;

    //- Null constructible
    constexpr null() noexcept {}

    //- Construct from Istream consumes no content.
    explicit constexpr null(Istream&) noexcept {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Global zero
static constexpr const zero Zero;

// IOstream Operators

//- Read from Istream consumes no content.
inline constexpr Istream& operator>>(Istream& is, zero&) noexcept
{
    return is;
}

//- Write to Ostream emits no content.
inline constexpr Ostream& operator<<(Ostream& os, const zero::null&) noexcept
{
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "zeroI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
