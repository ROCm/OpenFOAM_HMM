/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

InClass
    Foam

Description
    Spatial transformation functions for primitive fields.

SourceFiles
    transformList.C

\*---------------------------------------------------------------------------*/

#ifndef transformList_H
#define transformList_H

#include "transform.H"
#include "List.H"
#include "Map.H"
#include "EdgeMap.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Extend transform to work on list.
template<class T>
List<T> transform
(
    const tensor& rotTensor,
    const UList<T>& field
);

//- Apply transformation to list.
template<class T>
void transformList(const tensor&, UList<T>&);

template<class T>
void transformList(const tensor&, Map<T>&);

template<class T>
void transformList(const tensor&, EdgeMap<T>&);


template<>
inline void transformList(const tensor&, UList<bool>&)
{}
template<>
inline void transformList(const tensor&, Map<bool>&)
{}
template<>
inline void transformList(const tensor&, EdgeMap<bool>&)
{}


template<>
inline void transformList(const tensor&, UList<label>&)
{}
template<>
inline void transformList(const tensor&, Map<label>&)
{}
template<>
inline void transformList(const tensor&, EdgeMap<label>&)
{}


template<>
inline void transformList(const tensor&, UList<scalar>&)
{}
template<>
inline void transformList(const tensor&, Map<scalar>&)
{}
template<>
inline void transformList(const tensor&, EdgeMap<scalar>&)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "transformList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
