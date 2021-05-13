/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "PrimitivePatch.H"
#include "HashSet.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FaceList, class PointField>
Foam::labelList
Foam::PrimitivePatch<FaceList, PointField>::boundaryFaces() const
{
    // By definition boundary edges have a _single_ face attached,
    // but a face can easily have multiple boundary edges.

    const labelListList& edgeToFace = edgeFaces();

    labelHashSet bndFaces(2*nBoundaryEdges());

    for (label edgei = nInternalEdges(); edgei < edgeToFace.size(); ++edgei)
    {
        bndFaces.insert(edgeToFace[edgei][0]);
    }

    return bndFaces.sortedToc();
}


// ************************************************************************* //
