/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

InClass
    Foam::contiguous

Description
    Template function to specify if the data of a type are contiguous.

    The default function specifies that data are not contiguous.
    This is specialised for the types (eg, primitives) with contiguous data.

\*---------------------------------------------------------------------------*/

#ifndef contiguous_H
#define contiguous_H

#include "int.H"
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
template<class T, unsigned N> class FixedList;
template<class T> class Pair;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Default definition: (integral | floating-point) are contiguous
template<class T>
inline bool contiguous()
{
    return
    std::integral_constant
    <
        bool,
        std::is_integral<T>::value || std::is_floating_point<T>::value
    >::value;
}

//
// Fixed size containers of (integral | floating-point) are contiguous
//

template<>
inline bool contiguous<FixedList<bool, 2>>()               {return true;}
template<>
inline bool contiguous<Pair<bool>>()                       {return true;}

template<>
inline bool contiguous<FixedList<char, 2>>()               {return true;}
template<>
inline bool contiguous<Pair<char>>()                       {return true;}

template<>
inline bool contiguous<FixedList<int8_t, 2>>()             {return true;}
template<>
inline bool contiguous<Pair<int8_t>>()                     {return true;}

template<>
inline bool contiguous<FixedList<uint8_t, 2>>()            {return true;}
template<>
inline bool contiguous<Pair<uint8_t>>()                    {return true;}

template<>
inline bool contiguous<FixedList<int16_t, 2>>()            {return true;}
template<>
inline bool contiguous<Pair<int16_t>>()                    {return true;}

template<>
inline bool contiguous<FixedList<uint16_t, 2>>()           {return true;}
template<>
inline bool contiguous<Pair<uint16_t>>()                   {return true;}

template<>
inline bool contiguous<FixedList<int32_t, 2>>()            {return true;}
template<>
inline bool contiguous<Pair<int32_t>>()                    {return true;}

template<>
inline bool contiguous<FixedList<uint32_t, 2>>()           {return true;}
template<>
inline bool contiguous<Pair<uint32_t>>()                   {return true;}

template<>
inline bool contiguous<FixedList<int64_t, 2>>()            {return true;}
template<>
inline bool contiguous<Pair<int64_t>>()                    {return true;}

template<>
inline bool contiguous<FixedList<uint64_t, 2>>()           {return true;}
template<>
inline bool contiguous<Pair<uint64_t>>()                   {return true;}

template<>
inline bool contiguous<FixedList<float, 2>>()              {return true;}
template<>
inline bool contiguous<Pair<float>>()                      {return true;}

template<>
inline bool contiguous<FixedList<double, 2>>()             {return true;}
template<>
inline bool contiguous<Pair<double>>()                     {return true;}

template<>
inline bool contiguous<FixedList<long double, 2>>()        {return true;}
template<>
inline bool contiguous<Pair<long double>>()                {return true;}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
