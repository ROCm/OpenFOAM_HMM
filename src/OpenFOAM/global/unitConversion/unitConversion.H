/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

InNamespace
    Foam

Description
    Unit conversion functions

\*---------------------------------------------------------------------------*/

#ifndef unitConversion_H
#define unitConversion_H

#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Conversion from degrees to radians
inline constexpr scalar degToRad(const scalar deg) noexcept
{
    return (deg*M_PI/180.0);
}

//- Conversion from radians to degrees
inline constexpr scalar radToDeg(const scalar rad) noexcept
{
    return (rad*180.0/M_PI);
}

//- Multiplication factor for degrees to radians conversion
inline constexpr scalar degToRad() noexcept
{
    return (M_PI/180.0);
}

//- Multiplication factor for radians to degrees conversion
inline constexpr scalar radToDeg() noexcept
{
    return (180.0/M_PI);
}

//- Conversion from atm to Pa
inline constexpr scalar atmToPa(const scalar atm) noexcept
{
    return (atm*101325.0);
}

//- Conversion from atm to Pa
inline constexpr scalar paToAtm(const scalar pa) noexcept
{
    return (pa/101325.0);
}


//- User literal for degrees to radians conversion (integers)
inline constexpr scalar operator "" _deg(unsigned long long int deg) noexcept
{
    return (deg*M_PI/180.0);
}

//- User literal for degrees to radians conversion (floats)
inline constexpr scalar operator "" _deg(long double deg) noexcept
{
    return (deg*M_PI/180.0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
