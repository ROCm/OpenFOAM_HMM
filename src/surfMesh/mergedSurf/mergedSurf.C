/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2019 OpenCFD Ltd.
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

#include "mergedSurf.H"
#include "PatchTools.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mergedSurf::mergedSurf
(
    const meshedSurf& unmergedSurface,
    const scalar mergeDim
)
:
    mergedSurf()
{
    merge(unmergedSurface, mergeDim);
}


Foam::mergedSurf::mergedSurf
(
    const pointField& unmergedPoints,
    const faceList& unmergedFaces,
    const scalar mergeDim
)
:
    mergedSurf()
{
    merge(unmergedPoints, unmergedFaces, mergeDim);
}


Foam::mergedSurf::mergedSurf
(
    const pointField& unmergedPoints,
    const faceList& unmergedFaces,
    const labelList& originalIds,
    const scalar mergeDim
)
:
    mergedSurf()
{
    merge(unmergedPoints, unmergedFaces, originalIds, mergeDim);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::mergedSurf::use()
{
    return Pstream::parRun();
}


void Foam::mergedSurf::clear()
{
    points_.clear();
    faces_.clear();
    zones_.clear();
    pointsMap_.clear();
}


bool Foam::mergedSurf::merge
(
    const meshedSurf& unmergedSurface,
    const scalar mergeDim
)
{
    return
        merge
        (
            unmergedSurface.points(),
            unmergedSurface.faces(),
            unmergedSurface.zoneIds(),
            mergeDim
        );
}


bool Foam::mergedSurf::merge
(
    const pointField& unmergedPoints,
    const faceList& unmergedFaces,
    const scalar mergeDim
)
{
    return merge(unmergedPoints, unmergedFaces, labelList(), mergeDim);
}


bool Foam::mergedSurf::merge
(
    const pointField& unmergedPoints,
    const faceList& unmergedFaces,
    const labelList& originalIds,
    const scalar mergeDim
)
{
    if (!use())
    {
        clear();   // Extra safety?
        return false;
    }

    PatchTools::gatherAndMerge
    (
        mergeDim,
        primitivePatch
        (
            SubList<face>(unmergedFaces, unmergedFaces.size()),
            unmergedPoints
        ),
        points_,
        faces_,
        pointsMap_
    );


    // Now handle zone/region information

    globalIndex::gatherOp(originalIds, zones_);

    return true;
}


// ************************************************************************* //
