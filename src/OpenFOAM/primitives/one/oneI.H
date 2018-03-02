/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class arg2>
class innerProduct<one, arg2>
{
public:

    typedef arg2 type;
};


inline constexpr const one& operator*(const one& o, const one&) noexcept
{
    return o;
}

template<class Type>
inline constexpr const Type& operator*(const Type& val, const one&) noexcept
{
    return val;
}

template<class Type>
inline constexpr const Type& operator*(const one&, const Type& val) noexcept
{
    return val;
}

template<class Type>
inline constexpr const Type& operator&(const one&, const Type& val) noexcept
{
    return val;
}

inline constexpr const one& operator/(const one& o, const one&) noexcept
{
    return o;
}

template<class Type>
inline Type operator/(const one&, const Type& val)
{
    return scalar(1)/val;
}

template<class Type>
inline constexpr const Type& operator/(const Type& val, const one&) noexcept
{
    return val;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
