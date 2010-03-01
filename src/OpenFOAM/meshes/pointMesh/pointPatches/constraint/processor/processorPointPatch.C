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

#include "processorPointPatch.H"
#include "pointBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "pointMesh.H"
#include "faceList.H"
#include "primitiveFacePatch.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(processorPointPatch, 0);

addToRunTimeSelectionTable
(
    facePointPatch,
    processorPointPatch,
    polyPatch
);


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::processorPointPatch::initGeometry(PstreamBuffers& pBufs)
{
    // Algorithm:
    // Depending on whether the patch is a master or a slave, get the primitive
    // patch points and filter away the points from the global patch.

    // Create the reversed patch and pick up its points
    // so that the order is correct
    const polyPatch& pp = patch();

    faceList masterFaces(pp.size());

    forAll (pp, faceI)
    {
        masterFaces[faceI] = pp[faceI].reverseFace();
    }

    reverseMeshPoints_ = primitiveFacePatch
    (
        masterFaces,
        pp.points()
    ).meshPoints();
}


void Foam::processorPointPatch::calcGeometry(PstreamBuffers& pBufs)
{
    if (Pstream::parRun())
    {
        const boolList& collocated = procPolyPatch_.collocated();

        if (collocated.size() == 0)
        {
            separatedPoints_.setSize(0);
        }
        else if (collocated.size() == 1)
        {
            // Uniformly
            if (collocated[0])
            {
                separatedPoints_.setSize(0);
            }
            else
            {
                separatedPoints_ = identity(size());
            }
        }
        else
        {
            // Per face collocated or not.
            const labelListList& pointFaces = procPolyPatch_.pointFaces();

            DynamicList<label> separated;
            forAll(pointFaces, pfi)
            {
                if (!collocated[pointFaces[pfi][0]])
                {
                    separated.append(pfi);
                }
            }
            separatedPoints_.transfer(separated);
        }
    }

    if (debug)
    {
        Pout<< "processor:" << name()
            << " separated:" << separatedPoints_.size()
            << " out of points:" << size() << endl;
    }
}


void processorPointPatch::initMovePoints(PstreamBuffers&, const pointField&)
{}


void processorPointPatch::movePoints(PstreamBuffers&, const pointField&)
{}


void processorPointPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    facePointPatch::initUpdateMesh(pBufs);
    processorPointPatch::initGeometry(pBufs);
}


void processorPointPatch::updateMesh(PstreamBuffers& pBufs)
{
    facePointPatch::updateMesh(pBufs);
    processorPointPatch::calcGeometry(pBufs);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

processorPointPatch::processorPointPatch
(
    const polyPatch& patch,
    const pointBoundaryMesh& bm
)
:
    coupledFacePointPatch(patch, bm),
    procPolyPatch_(refCast<const processorPolyPatch>(patch))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

processorPointPatch::~processorPointPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const labelList& processorPointPatch::reverseMeshPoints() const
{
    return reverseMeshPoints_;
}


const labelList& processorPointPatch::separatedPoints() const
{
    return separatedPoints_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
