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

#include "mergedSurf.H"
#include "PatchTools.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mergedSurf::mergedSurf()
:
    points_(),
    faces_(),
    zones_(),
    pointsMap_()
{}


Foam::mergedSurf::mergedSurf
(
    const meshedSurf& surf,
    const scalar mergeDim
)
:
    mergedSurf()
{
    merge(surf, mergeDim);
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
    const meshedSurf& surf,
    const scalar mergeDim
)
{
    // needed for extra safety?
    // clear();

    if (!use())
    {
        return false;
    }

    PatchTools::gatherAndMerge
    (
        mergeDim,
        primitivePatch
        (
            SubList<face>(surf.faces(), surf.faces().size()),
            surf.points()
        ),
        points_,
        faces_,
        pointsMap_
    );

    // Now handle zone/region information
    List<labelList> allZones(Pstream::nProcs());
    allZones[Pstream::myProcNo()] = surf.zoneIds();
    Pstream::gatherList(allZones);

    if (Pstream::master())
    {
        zones_ = ListListOps::combine<labelList>
        (
            allZones,
            accessOp<labelList>()
        );
    }
    allZones.clear();

    return true;
}


// ************************************************************************* //
