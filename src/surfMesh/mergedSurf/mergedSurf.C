/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    const labelList& origZoneIds,
    const labelList& origFaceIds,
    const scalar mergeDim
)
:
    mergedSurf()
{
    merge
    (
        unmergedPoints,
        unmergedFaces,
        origZoneIds,
        origFaceIds,
        mergeDim
    );
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
    pointsMap_.clear();

    zoneIds_.clear();
    faceIds_.clear();
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
            unmergedSurface.faceIds(),
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
    return
        merge
        (
            unmergedPoints,
            unmergedFaces,
            labelList(),
            labelList(),
            mergeDim
        );
}


bool Foam::mergedSurf::merge
(
    const pointField& unmergedPoints,
    const faceList& unmergedFaces,
    const labelList& origZoneIds,
    const labelList& origFaceIds,
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


    // Now handle per-face information

    globalIndex::gatherOp(origZoneIds, zoneIds_);
    globalIndex::gatherOp(origFaceIds, faceIds_);

    return true;
}


// ************************************************************************* //
