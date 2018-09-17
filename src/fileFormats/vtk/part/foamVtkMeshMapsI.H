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

\*---------------------------------------------------------------------------*/

#include "foamVtkMeshMaps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::foamVtkMeshMaps::foamVtkMeshMaps(const label size)
:
    cellMap_(size),
    pointMap_(size),
    additionalIds_(size)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

inline Foam::foamVtkMeshMaps::~foamVtkMeshMaps()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::foamVtkMeshMaps::clear()
{
    cellMap_.clear();
    pointMap_.clear();
    additionalIds_.clear();
}


inline const Foam::labelList&
Foam::foamVtkMeshMaps::cellMap() const
{
    return cellMap_;
}


inline Foam::DynamicList<Foam::label>&
Foam::foamVtkMeshMaps::cellMap()
{
    return cellMap_;
}


inline const Foam::labelList&
Foam::foamVtkMeshMaps::pointMap() const
{
    return pointMap_;
}


inline Foam::DynamicList<Foam::label>&
Foam::foamVtkMeshMaps::pointMap()
{
    return pointMap_;
}


inline const Foam::labelList&
Foam::foamVtkMeshMaps::additionalIds() const
{
    return additionalIds_;
}


inline Foam::DynamicList<Foam::label>&
Foam::foamVtkMeshMaps::additionalIds()
{
    return additionalIds_;
}


// ************************************************************************* //
