/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include <cstdio>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// NOTE: with C++11 could consider variadic templates for a more general
// sprintf implementation

template<class PrimitiveType>
Foam::word Foam::stringOps::name
(
    const char* fmt,
    const PrimitiveType& val
)
{
    // same concept as GNU/BSD asprintf()
    // use snprintf with zero to determine the number of characters required

    int n = ::snprintf(0, 0, fmt, val);
    if (n > 0)
    {
        char buf[n+1];
        ::snprintf(buf, n+1, fmt, val);
        buf[n] = 0;

        // no stripping desired
        return word(buf, false);
    }

    return word::null;
}


template<class PrimitiveType>
Foam::word Foam::stringOps::name
(
    const std::string& fmt,
    const PrimitiveType& val
)
{
    return stringOps::name(fmt.c_str(), val);
}


// ************************************************************************* //
