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

\*---------------------------------------------------------------------------*/

#include "cyclicPointPatch.H"
#include "pointBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "pointMesh.H"
#include "edgeList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(cyclicPointPatch, 0);

addToRunTimeSelectionTable
(
    facePointPatch,
    cyclicPointPatch,
    polyPatch
);


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::cyclicPointPatch::initGeometry(PstreamBuffers&)
{}


void Foam::cyclicPointPatch::calcGeometry(PstreamBuffers&)
{
    const edgeList& cp = cyclicPolyPatch_.coupledPoints();
    const labelList& mp = cyclicPolyPatch_.meshPoints();

    DynamicList<label> separated;
    forAll(cp, i)
    {
        const edge& coupledSet = cp[i];

        // Assume all points are separated.
        separated.append(coupledSet[0]);
        separated.append(coupledSet[1]);
    }
    separatedPoints_.transfer(separated);

    if (debug)
    {
        Pout<< "cyclic:" << cyclicPolyPatch_.name()
            << " separated:" << separatedPoints_.size()
            << " out of points:" << mp.size() << endl;
    }
}


void cyclicPointPatch::initMovePoints(PstreamBuffers&, const pointField&)
{}


void cyclicPointPatch::movePoints(PstreamBuffers&, const pointField&)
{}


void cyclicPointPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    facePointPatch::initUpdateMesh(pBufs);
    cyclicPointPatch::initGeometry(pBufs);
}


void cyclicPointPatch::updateMesh(PstreamBuffers& pBufs)
{
    facePointPatch::updateMesh(pBufs);
    cyclicPointPatch::calcGeometry(pBufs);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

cyclicPointPatch::cyclicPointPatch
(
    const polyPatch& patch,
    const pointBoundaryMesh& bm
)
:
    coupledFacePointPatch(patch, bm),
    cyclicPolyPatch_(refCast<const cyclicPolyPatch>(patch))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cyclicPointPatch::~cyclicPointPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const edgeList& cyclicPointPatch::transformPairs() const
{
    return cyclicPolyPatch_.coupledPoints();
}


const labelList& cyclicPointPatch::separatedPoints() const
{
    return separatedPoints_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
