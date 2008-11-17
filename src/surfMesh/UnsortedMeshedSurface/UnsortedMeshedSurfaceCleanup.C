/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "UnsortedMeshedSurface.H"
#include "mergePoints.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Remove badly degenerate faces, double faces.
template<class Face>
void Foam::UnsortedMeshedSurface<Face>::cleanup(const bool verbose)
{
    // merge points (already done for STL, TRI)
    stitchFaces(SMALL, verbose);

    checkFaces(verbose);
    ParentType::checkEdges(verbose);
}


template<class Face>
bool Foam::UnsortedMeshedSurface<Face>::stitchFaces
(
    const scalar tol,
    const bool verbose
)
{
    List<label> faceMap;
    bool hasMerged = ParentType::stitchFaces(faceMap, tol, verbose);

    remapRegions(faceMap);
    return hasMerged;
}


// Remove badly degenerate faces and double faces.
template<class Face>
bool Foam::UnsortedMeshedSurface<Face>::checkFaces(const bool verbose)
{
    List<label> faceMap;
    bool changed = ParentType::checkFaces(faceMap, verbose);

    remapRegions(faceMap);
    return changed;
}


template<class Face>
Foam::label Foam::UnsortedMeshedSurface<Face>::triangulate()
{
    List<label> faceMap;
    label nTri = ParentType::triangulate(this->storedFaces(), faceMap);

    remapRegions(faceMap);
    return nTri;
}

// ************************************************************************* //
