/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline constexpr zero operator+(const zero&, const zero&) noexcept
{
    return Zero;
}

template<class Type>
inline constexpr const Type& operator+(const Type& val, const zero&) noexcept
{
    return val;
}

template<class Type>
inline constexpr const Type& operator+(const zero&, const Type& val) noexcept
{
    return val;
}

inline constexpr zero operator-(const zero&, const zero&) noexcept
{
    return Zero;
}

template<class Type>
inline constexpr const Type& operator-(const Type& val, const zero&) noexcept
{
    return val;
}

template<class Type>
inline constexpr Type operator-(const zero&, const Type& val)
{
    return -val;
}

inline constexpr zero operator*(const zero&, const zero&) noexcept
{
    return Zero;
}

template<class Type>
inline constexpr zero operator*(const Type& val, const zero&) noexcept
{
    return Zero;
}

template<class Type>
inline constexpr zero operator*(const zero&, const Type& val) noexcept
{
    return Zero;
}

template<class Type>
inline constexpr zero operator/(const zero&, const Type& val) noexcept
{
    return Zero;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
