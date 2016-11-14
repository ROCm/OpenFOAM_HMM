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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// A triFace surface only handles triangulated faces
template<>
inline bool MeshedSurface<triFace>::isTri()
{
    return true;
}


// A labelledTri surface only handles triangulated faces
template<>
inline bool MeshedSurface<labelledTri>::isTri()
{
    return true;
}


// Number of triangles for a triFace surface
template<>
inline label MeshedSurface<triFace>::nTriangles() const
{
    return ParentType::size();
}

// Number of triangles for a labelledTri surface
template<>
inline label MeshedSurface<labelledTri>::nTriangles() const
{
    return ParentType::size();
}


// Inplace triangulation of triFace surface = no-op
template<>
inline label MeshedSurface<triFace>::triangulate()
{
    return 0;
}

// Inplace triangulation of labelledTri surface = no-op
template<>
inline label MeshedSurface<labelledTri>::triangulate()
{
    return 0;
}

// Inplace triangulation of triFace surface (with face map) = no-op
template<>
inline label MeshedSurface<triFace>::triangulate(List<label>& faceMap)
{
    if (notNull(faceMap))
    {
        faceMap.clear();
    }
    return 0;
}

// Inplace triangulation of labelledTri surface (with face map) = no-op
template<>
inline label MeshedSurface<labelledTri>::triangulate(List<label>& faceMap)
{
    if (notNull(faceMap))
    {
        faceMap.clear();
    }
    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
