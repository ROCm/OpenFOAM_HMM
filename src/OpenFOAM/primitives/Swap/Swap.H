/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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
    Swap arguments as per std::swap, but in Foam namespace.

\*---------------------------------------------------------------------------*/

#ifndef Swap_H
#define Swap_H

#include <type_traits>
#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Swap non-array types as per std::swap example, but in Foam namespace.
//  \sa http://www.cplusplus.com/reference/utility/swap/
template<class T>
void Swap(T& a, T& b)
{
    // compile-time resolution with std::enable_if not yet working
    if (std::is_fundamental<T>::value || std::is_pointer<T>::value)
    {
        // Use copy/assign for simple types
        const T tmp = a;
        a = b;
        b = tmp;
    }
    else
    {
        // Use move/assignment
        T tmp(std::move(a));
        a = std::move(b);
        b = std::move(tmp);
    }
}


//- Swap array types as per std::swap example, but in Foam namespace.
//  \sa http://www.cplusplus.com/reference/utility/swap/
template<class T, size_t N>
void Swap(T (&a)[N], T (&b)[N])
{
    for (size_t i = 0; i < N; ++i)
    {
        Foam::Swap(a[i], b[i]);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
