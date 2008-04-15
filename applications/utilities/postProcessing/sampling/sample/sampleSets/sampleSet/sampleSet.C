/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "sampleSet.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "meshSearch.H"
#include "writer.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

const scalar sampleSet::tol = 1e-6;

defineTypeNameAndDebug(sampleSet, 0);
defineRunTimeSelectionTable(sampleSet, word);

autoPtr<sampleSet> sampleSet::New
(
    const word& sampleType,
    const polyMesh& mesh,
    meshSearch& searchEngine,
    const dictionary& dict
)
{
    wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_
            ->find(sampleType);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "sampleSet::New(const word&, "
            "const polyMesh&, meshSearch&, const dictionary&)"
        )   << "Unknown sample type " << sampleType
            << endl << endl
            << "Valid sample types : " << endl
            << wordConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<sampleSet>
    (
        cstrIter()
        (
            mesh,
            searchEngine,
            dict
        )
    );
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::sampleSet::getBoundaryCell(const label faceI) const
{
    return mesh().faceOwner()[faceI];
}


Foam::label Foam::sampleSet::getCell
(
    const label faceI,
    const point& sample
) const
{
    if (faceI == -1)
    {
        FatalErrorIn
        (
            "sampleSet::getCell(const label, const point&)"
        )   << "Illegal face label " << faceI
            << abort(FatalError);
    }

    if (faceI >= mesh().nInternalFaces())
    {
        label cellI = getBoundaryCell(faceI);

        if (!mesh().pointInCell(sample, cellI))
        {
            FatalErrorIn
            (
                "sampleSet::getCell(const label, const point&)"
            )   << "Found cell " << cellI << " using face " << faceI
                << ". But cell does not contain point " << sample
                << abort(FatalError);
        }
        return cellI;
    }
    else
    {
        // Try owner and neighbour to see which one contains sample

        label cellI = mesh().faceOwner()[faceI];

        if (mesh().pointInCell(sample, cellI))
        {
            return cellI;
        }
        else
        {
            cellI = mesh().faceNeighbour()[faceI];

            if (mesh().pointInCell(sample, cellI))
            {
                return cellI;
            }
            else
            {
                FatalErrorIn
                (
                    "sampleSet::getCell(const label, const point&)"
                )   << "None of the neighbours of face "
                    << faceI << " contains point " << sample
                    << abort(FatalError);

                return -1;
            }
        }
    }
}


Foam::scalar Foam::sampleSet::calcSign
(
    const label faceI,
    const point& sample
) const
{
    vector vec = sample - mesh().faceCentres()[faceI];

    scalar magVec = mag(vec);

    if (magVec < VSMALL)
    {
        // sample on face centre. Regard as inside
        return -1;
    }

    vec /= magVec;

    vector n = mesh().faceAreas()[faceI];

    n /= mag(n) + VSMALL;

    return n & vec;    
}


// Return face (or -1) of face which is within smallDist of sample
Foam::label Foam::sampleSet::findNearFace
(
    const label cellI,
    const point& sample,
    const scalar smallDist
) const
{
    const cell& myFaces = mesh().cells()[cellI];

    forAll(myFaces, myFaceI)
    {
        const face& f = mesh().faces()[myFaces[myFaceI]];

        pointHit inter = f.nearestPoint(sample, mesh().points());

        scalar dist;

        if (inter.hit())
        {
            dist = mag(inter.hitPoint() - sample);
        }
        else
        {
            dist = mag(inter.missPoint() - sample);
        }

        if (dist < smallDist)
        {
            return myFaces[myFaceI];
        }
    }
    return -1;
}


// 'Pushes' point facePt (which is almost on face) in direction of cell centre
// so it is clearly inside.
Foam::point Foam::sampleSet::pushIn
(
    const point& facePt,
    const label faceI
) const
{
    label cellI = mesh().faceOwner()[faceI];

    const point& cellCtr = mesh().cellCentres()[cellI];

    point newSample =
        facePt + tol*(cellCtr - facePt);

    if (!searchEngine().pointInCell(newSample, cellI))
    {
        FatalErrorIn
        (
            "sampleSet::pushIn(const point&, const label)"
        )   << "After pushing " << facePt << " to " << newSample
            << " it is still outside faceI " << faceI << endl
            << "Please change your starting point"
            << abort(FatalError);
    }
    //Info<< "pushIn : moved " << facePt << " to " << newSample
    //    << endl;

    return newSample;
}


// Calculates start of tracking given samplePt and first boundary intersection
// (bPoint, bFaceI). bFaceI == -1 if no boundary intersection.
// Returns true if trackPt is sampling point
bool Foam::sampleSet::getTrackingPoint
(
    const vector& offset,
    const point& samplePt,
    const point& bPoint,
    const label bFaceI,

    point& trackPt,
    label& trackCellI,
    label& trackFaceI
) const
{
    const scalar smallDist = mag(tol*offset);

    bool isGoodSample = false;

    if (bFaceI == -1)
    {
        // No boundary intersection. Try and find cell samplePt is in
        trackCellI = mesh().findCell(samplePt);

        if
        (
            (trackCellI == -1)
         || !mesh().pointInCell(samplePt, trackCellI)
        )
        {
            // Line samplePt - end_ does not intersect domain at all.
            // (or is along edge)
            //Info<< "getTrackingPoint : samplePt outside domain : "
            //    << "  samplePt:" << samplePt
            //    << endl;

            trackCellI = -1;
            trackFaceI = -1;

            isGoodSample = false;
        }
        else
        {
            // start is inside. Use it as tracking point
            //Info<< "getTrackingPoint : samplePt inside :"
            //    << "  samplePt:" << samplePt
            //    << "  trackCellI:" << trackCellI
            //    << endl;

            trackPt = samplePt;
            trackFaceI = -1;

            isGoodSample = true;
        }
    }
    else if (mag(samplePt - bPoint) < smallDist)
    {
        //Info<< "getTrackingPoint : samplePt:" << samplePt
        //    << " close to bPoint:"
        //    << bPoint << endl;

        // samplePt close to bPoint. Snap to it
        trackPt = pushIn(bPoint, bFaceI);
        trackFaceI = bFaceI;
        trackCellI = getBoundaryCell(trackFaceI);

        isGoodSample = true;
    }
    else
    {
        scalar sign = calcSign(bFaceI, samplePt);

        if (sign < 0)
        {
            // samplePt inside or marginally outside.
            trackPt = samplePt;
            trackFaceI = -1;
            trackCellI = mesh().findCell(trackPt);

            isGoodSample = true;
        }
        else
        {
            // samplePt outside. use bPoint
            trackPt = pushIn(bPoint, bFaceI);
            trackFaceI = bFaceI;
            trackCellI = getBoundaryCell(trackFaceI);

            isGoodSample = false;
        }
    }

    if (debug)
    {
        Info<< "sampleSet::getTrackingPoint :"
            << " offset:" << offset
            << " samplePt:" << samplePt
            << " bPoint:" << bPoint
            << " bFaceI:" << bFaceI
            << endl << "   Calculated first tracking point :"
            << " trackPt:" << trackPt
            << " trackCellI:" << trackCellI
            << " trackFaceI:" << trackFaceI
            << " isGoodSample:" << isGoodSample
            << endl;
    }

    return isGoodSample;
}            


void Foam::sampleSet::setSamples
(
    const List<point>& samplingPts,
    const labelList& samplingCells,
    const labelList& samplingFaces,
    const labelList& samplingSegments,
    const scalarList& samplingCurveDist
)
{
    setSize(samplingPts.size());
    cells_.setSize(samplingCells.size());
    faces_.setSize(samplingFaces.size());
    segments_.setSize(samplingSegments.size());
    curveDist_.setSize(samplingCurveDist.size());

    if
    (
        (cells_.size() != size())
     || (faces_.size() != size())
     || (segments_.size() != size())
     || (curveDist_.size() != size())
    )
    {
        FatalErrorIn("sampleSet::setSamples()")
            << "sizes not equal : "
            << "  points:" << size()
            << "  cells:" << cells_.size()
            << "  faces:" << faces_.size()
            << "  segments:" << segments_.size()
            << "  curveDist:" << curveDist_.size()
            << abort(FatalError);
    }

    forAll(samplingPts, sampleI)
    {
        operator[](sampleI) = samplingPts[sampleI];
    }
    cells_ = samplingCells;
    faces_ = samplingFaces;
    segments_ = samplingSegments;
    curveDist_ = samplingCurveDist;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh, name
Foam::sampleSet::sampleSet
(
    const polyMesh& mesh,
    meshSearch& searchEngine,
    const word& name,
    const word& axis
)
:
    coordSet(name, axis),
    mesh_(mesh),
    searchEngine_(searchEngine),
    segments_(0),
    curveDist_(0),
    cells_(0),
    faces_(0)
{}


// Construct from dictionary
Foam::sampleSet::sampleSet
(
    const polyMesh& mesh,
    meshSearch& searchEngine,
    const dictionary& dict          
)
:
    coordSet(dict.lookup("name"), dict.lookup("axis")),
    mesh_(mesh),
    searchEngine_(searchEngine),
    segments_(0),
    curveDist_(0),
    cells_(0),
    faces_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampleSet::~sampleSet()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Ostream& Foam::sampleSet::write(Ostream& os) const
{
    coordSet::write(os);

    os  << endl << "\t(cellI)\t(faceI)"
        << endl;
    forAll(*this, sampleI)
    {
        os  << '\t' << cells_[sampleI]
            << '\t' << faces_[sampleI]
            << endl;
    }
    return os;
}

// ************************************************************************* //
