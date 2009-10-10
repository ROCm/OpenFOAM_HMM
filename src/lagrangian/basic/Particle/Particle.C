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

#include "Particle.H"
#include "Cloud.H"
#include "wedgePolyPatch.H"
#include "symmetryPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "processorPolyPatch.H"
#include "wallPolyPatch.H"
#include "transform.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ParticleType>
void Foam::Particle<ParticleType>::findFaces
(
    const vector& endPosition,
    DynamicList<label>& faceList
) const
{
    const polyMesh& mesh = cloud_.polyMesh_;
    const labelList& faces = mesh.cells()[celli_];
    const vector& C = mesh.cellCentres()[celli_];

    faceList.clear();

    forAll(faces, i)
    {
        label facei = faces[i];
        scalar lam = lambda(C, endPosition, facei);

        if ((lam > 0) && (lam < 1.0))
        {
            faceList.append(facei);
        }
    }
}


template<class ParticleType>
void Foam::Particle<ParticleType>::findFaces
(
    const vector& endPosition,
    const label celli,
    const scalar stepFraction,
    DynamicList<label>& faceList
) const
{
    const polyMesh& mesh = cloud_.pMesh();
    const labelList& faces = mesh.cells()[celli];
    const vector& C = mesh.cellCentres()[celli];

    faceList.clear();
    forAll(faces, i)
    {
        label facei = faces[i];
        scalar lam = lambda(C, endPosition, facei, stepFraction);

        if ((lam > 0) && (lam < 1.0))
        {
            faceList.append(facei);
        }
    }
}


template<class ParticleType>
bool Foam::Particle<ParticleType>::insideCellExact
(
    const vector& testPt,
    const label celli,
    bool beingOnAFaceMeansOutside
) const
{
    const polyMesh& mesh = cloud_.pMesh();
    const labelList& faces = mesh.cells()[celli];
    const vector& C = mesh.cellCentres()[celli];

    label nFaceCrossings = 0;

    // The vector from the cell centre to the end point
    vector delta = testPt - C;

    forAll (faces, i)
    {
        label facei = faces[i];

        pointHit inter = mesh.faces()[facei].intersection
        (
            C,
            delta,
            mesh.faceCentres()[facei],
            mesh.points(),
            intersection::HALF_RAY,
            Cloud<ParticleType>::intersectionTolerance
        );

        if (inter.hit())
        {
            // Pout<< "insideCellExact cell " << celli
            //     << " face " << facei << " "
            //     << inter.distance() << endl;

            if (beingOnAFaceMeansOutside)
            {
                if (inter.distance() <= 1.0)
                {
                    // This face was actually crossed.

                    nFaceCrossings++;
                }
            }
            else
            {
                if (inter.distance() < 1.0)
                {
                    // This face was actually crossed.

                    nFaceCrossings++;
                }
            }
        }
    }

    if (nFaceCrossings > 1)
    {
        Pout<< "In cell " << celli_ << " there were " << nFaceCrossings
            << " face crossings detected tracking from concave cell centre to "
            << " endPosition"
            << endl;
    }

    if (nFaceCrossings % 2 == 0)
    {
        // Even number of face crossings, so the testPt must be in the
        // cell.

        return true;
    }

    return false;
}


template<class ParticleType>
template<class TrackData>
void Foam::Particle<ParticleType>::trackToFaceConcave
(
    scalar& trackFraction,
    const vector& endPosition,
    TrackData& td
)
{
    facei_ = -1;

    const polyMesh& mesh = cloud_.pMesh();
    const labelList& faces = mesh.cells()[celli_];

    // Check all possible face crossings to see if they are actually
    // crossed, determining if endPosition is outside the current
    // cell.  This allows situations where the cell is outside the
    // cell to start with and enters the cell at the end of the track
    // to be identified.

    // Pout<< nl << "Outside test:" << endl;

    if (insideCellExact(endPosition, celli_, false))
    {
        // Even number of face crossings, so the particle must end up
        // still in the cell.

        position_ = endPosition;
        trackFraction = 1.0;
        return;
    }

    // Pout<< nl << origProc_ << " "
    //     << origId_ << " "
    //     << position_ << " "
    //     << endPosition << " "
    //     << stepFraction_ << " "
    //     << celli_
    //     << endl;

    // The particle *must* have left the cell.

    // a) It may have crossed a face not yet identified by testing
    //    faces using the cell centre to endPosition line, so the
    //    potentially crossed faces of the position to endPosition
    //    line must be assessed.

    // b) It may have been outside the cell in the first place, and, despite
    //    trying to pick up more faces using a) the correct face to be crossed
    //    is not knowable.  A best guess will be used, with the expectation that
    //    the tracking in the destination cell will be able to recover form a
    //    bad guess.

    // For all face assessments, a full intersection test is required,
    // as nothing can be assumed about the order of crossing the
    // planes of faces.

    const vector deltaPosition = endPosition - position_;
    vector deltaTrack =
        mag(mesh.cellCentres()[celli_] - position_)
       *deltaPosition/(mag(deltaPosition) + VSMALL);

    // Pout<< "Inside test:" << endl;

    if (insideCellExact(position_, celli_, false))
    {
        // Pout<< "The particle starts inside the cell and ends up outside of it"
        //     << nl << position_ << " " << position_ + deltaTrack
        //     << endl;

        // The particle started inside the cell and finished outside
        // of it, find which face to cross

        scalar tmpLambda = GREAT;
        scalar correctLambda = GREAT;

        forAll(faces, i)
        {
            label facei = faces[i];

            // Use exact intersection.

            // TODO: A correction is required for moving meshes to
            // calculate the correct lambda value.

            pointHit inter = mesh.faces()[facei].intersection
            (
                position_,
                deltaPosition,
                mesh.faceCentres()[facei],
                mesh.points(),
                intersection::HALF_RAY,
                Cloud<ParticleType>::intersectionTolerance
            );

            if (inter.hit())
            {
                tmpLambda = inter.distance();

                // Pout<< facei << " " << tmpLambda << endl;

                if
                (
                    tmpLambda <= 1.0
                 && tmpLambda < correctLambda
                )
                {
                    // This face is crossed before any other that has
                    // been found so far

                    correctLambda = tmpLambda;
                    facei_ = facei;
                }
            }
        }

        if (facei_ > -1)
        {
            if (!cloud_.internalFace(facei_))
            {
                // For a patch face, allow a small value of lambda to
                // ensure patch interactions occur.

                label patchi = patch(facei_);
                const polyPatch& patch = mesh.boundaryMesh()[patchi];

                if (isA<wallPolyPatch>(patch))
                {
                    if ((mesh.faceAreas()[facei_] & deltaPosition) <= 0)
                    {
                        // The particle has hit a wall face but it is
                        // heading in the wrong direction with respect to
                        // the face normal

                        // Do not trigger a face hit and move the position
                        // towards the cell centre

                        // Pout<< "Hit a wall face heading the wrong way"
                        //     << endl;

                        const point& cc = mesh.cellCentres()[celli_];
                        position_ +=
                            Cloud<ParticleType>::trackingRescueTolerance
                           *(cc - position_);

                        facei_ = -1;
                    }
                }
            }
            else
            {
                if (correctLambda < Cloud<ParticleType>::minValidTrackFraction)
                {
                    // The particle is not far enough away from the face
                    // to decide if it is valid crossing. Let it move a
                    // little without crossing the face to resolve the
                    // ambiguity.

                    // Pout<< "Ambiguous face crossing, correcting towards cell "
                    //     << "centre and not crossing face" << endl;

                    // const point& cc = mesh.cellCentres()[celli_];
                    // position_ +=
                    //     Cloud<ParticleType>::trackingRescueTolerance
                    //    *(cc - position_);

                    // Pout<< "Ambiguous face crossing. " << endl;

                    facei_ = -1;
                }

                // If the face hit was not on a wall, add a small
                // amount to the track to move it off the face, If it
                // was not an ambiguous face crossing, this makes sure
                // the face is not ambiguous next tracking step.  If
                // it was ambiguous, this should resolve it.

                correctLambda += Cloud<ParticleType>::minValidTrackFraction;
            }

            trackFraction = correctLambda;
            position_ += trackFraction*(endPosition - position_);
        }
        else
        {
            // Pout<< "Particle " << origProc_ << " " << origId_
            //     << " started inside cell " << celli_ << " and finished outside"
            //     << " of it, but did not find a face to cross"
            //     << endl;

            const point& cc = mesh.cellCentres()[celli_];
            position_ +=
                Cloud<ParticleType>::trackingRescueTolerance*(cc - position_);
        }
    }
    else
    {
        // Pout<< "The particle started outside of the cell" << endl;

        // Find which cell the particle should be in.

        const labelList& cPts = mesh.cellPoints(celli_);

        DynamicList<label> checkedCells;

        bool found = false;

        forAll(cPts, cPtI)
        {
            label ptI = cPts[cPtI];

            const labelList& pCs = mesh.pointCells(ptI);

            forAll(pCs, pCI)
            {
                label cellI = pCs[pCI];

                if (findIndex(checkedCells, cellI) == -1)
                {
                    checkedCells.append(cellI);

                    if (insideCellExact(position_, cellI, false))
                    {
                        found = true;

                        celli_ = cellI;

                        break;
                    }
                }
            }

            if (found)
            {
                break;
            }
        }

        if (!found)
        {
            // Pout<< "Didn't find a new cell after searching "
            //     << checkedCells << endl;

            const point& cc = mesh.cellCentres()[celli_];
            position_ +=
                Cloud<ParticleType>::trackingRescueTolerance*(cc - position_);
        }
        // else
        // {
        //     Pout<< "Found new cell " << celli_
        //         << " by searching " << checkedCells
        //         << endl;
        // }
    }

    // Pout<< facei_ << " " << celli_ << endl;

    if (facei_ > -1)
    {
        faceAction(trackFraction, endPosition, td);
    }

    // Pout<< facei_ << " " << celli_ << endl;
}


template<class ParticleType>
template<class TrackData>
void Foam::Particle<ParticleType>::trackToFaceConvex
(
    scalar& trackFraction,
    const vector& endPosition,
    TrackData& td
)
{
    facei_ = -1;

    DynamicList<label>& faces = cloud_.labels_;
    findFaces(endPosition, faces);

    if (faces.empty())
    {
        // endPosition is inside the cell

        position_ = endPosition;
        trackFraction = 1.0;
        return;
    }

    // A face has been hit

    scalar lambdaMin = GREAT;

    if (faces.size() == 1)
    {
        lambdaMin = lambda(position_, endPosition, faces[0], stepFraction_);
        facei_ = faces[0];
    }
    else
    {
        forAll(faces, i)
        {
            scalar lam =
                lambda(position_, endPosition, faces[i], stepFraction_);

            if (lam < lambdaMin)
            {
                lambdaMin = lam;
                facei_ = faces[i];
            }
        }
    }

    // For warped faces the particle can be 'outside' the cell.
    // This will yield a lambda larger than 1, or smaller than 0.

    // For values < 0, the particle travels away from the cell and we
    // don't move the particle (except by a small value to move it off
    // the face), only change cell.

    if (static_cast<ParticleType&>(*this).softImpact())
    {
        // Soft-sphere particles can travel outside the domain
        // but we don't use lambda since this the particle
        // is going away from face

        trackFraction = 1.0;
        position_ = endPosition;
    }
    else if (lambdaMin <= 0.0)
    {
        // Pout<< "convex tracking recovery "
        //     << origId_ << " "
        //     << origProc_ << " "
        //     << position_ << " "
        //     << endPosition << " "
        //     << stepFraction_ << " "
        //     << lambdaMin << " "
        //     << celli_ << " "
        //     << facei_ << " "
        //     << endl;

        trackFraction = Cloud<ParticleType>::trackingRescueTolerance;
        position_ += trackFraction*(endPosition - position_);
    }
    else
    {
        if (lambdaMin <= 1.0)
        {
            trackFraction = lambdaMin;
            position_ += trackFraction*(endPosition - position_);
        }
        else
        {
            // For values larger than 1, we move the particle to endPosition
            // only.
            trackFraction = 1.0;
            position_ = endPosition;
        }
    }

    faceAction(trackFraction, endPosition, td);

    // If the trackFraction = 0 something went wrong.
    // Either the particle is flipping back and forth across a face perhaps
    // due to velocity interpolation errors or it is in a "hole" in the mesh
    // caused by face warpage.
    // In both cases resolve the positional ambiguity by moving the particle
    // slightly towards the cell-centre.

    if (trackFraction < Cloud<ParticleType>::minValidTrackFraction)
    {
        // Pout<< "convex tracking error "
        //     << origId_ << " "
        //     << origProc_ << " "
        //     << position_ << " "
        //     << endPosition << " "
        //     << trackFraction << " "
        //     << stepFraction_ << " "
        //     << lambdaMin << " "
        //     << celli_ << " "
        //     << facei_ << " "
        //     << endl;

        const polyMesh& mesh = cloud_.pMesh();

        const point& cc = mesh.cellCentres()[celli_];
        position_ +=
            Cloud<ParticleType>::trackingRescueTolerance*(cc - position_);
    }
}

template<class ParticleType>
template<class TrackData>
void Foam::Particle<ParticleType>::faceAction
(
    scalar& trackFraction,
    const vector& endPosition,
    TrackData& td
)
{
    const polyMesh& mesh = cloud_.pMesh();

    if (cloud_.internalFace(facei_))
    {
        // Internal face, change cell

        if (celli_ == mesh.faceOwner()[facei_])
        {
            celli_ = mesh.faceNeighbour()[facei_];
        }
        else if (celli_ == mesh.faceNeighbour()[facei_])
        {
            celli_ = mesh.faceOwner()[facei_];
        }
        else
        {
            FatalErrorIn("Particle::trackToFace(const vector&, TrackData&)")
            << "addressing failure" << nl
                << abort(FatalError);
        }
    }
    else
    {
        ParticleType& p = static_cast<ParticleType&>(*this);

        // Soft-sphere algorithm ignores the boundary
        if (p.softImpact())
        {
            trackFraction = 1.0;
            position_ = endPosition;
        }

        label patchi = patch(facei_);
        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        if (!p.hitPatch(patch, td, patchi))
        {
            if (isA<wedgePolyPatch>(patch))
            {
                p.hitWedgePatch
                (
                    static_cast<const wedgePolyPatch&>(patch), td
                );
            }
            else if (isA<symmetryPolyPatch>(patch))
            {
                p.hitSymmetryPatch
                (
                    static_cast<const symmetryPolyPatch&>(patch), td
                );
            }
            else if (isA<cyclicPolyPatch>(patch))
            {
                p.hitCyclicPatch
                (
                    static_cast<const cyclicPolyPatch&>(patch), td
                );
            }
            else if (isA<processorPolyPatch>(patch))
            {
                p.hitProcessorPatch
                (
                    static_cast<const processorPolyPatch&>(patch), td
                );
            }
            else if (isA<wallPolyPatch>(patch))
            {
                p.hitWallPatch
                (
                    static_cast<const wallPolyPatch&>(patch), td
                );
            }
            else
            {
                p.hitPatch(patch, td);
            }
        }
    }
}

template<class ParticleType>
template<class TrackData>
void Foam::Particle<ParticleType>::prepareForParallelTransfer
(
    const label patchi,
    TrackData& td
)
{
    // Convert the face index to be local to the processor patch
    facei_ = patchFace(patchi, facei_);
}


template<class ParticleType>
template<class TrackData>
void Foam::Particle<ParticleType>::correctAfterParallelTransfer
(
    const label patchi,
    TrackData& td
)
{
    const processorPolyPatch& ppp =
        refCast<const processorPolyPatch>
        (cloud_.pMesh().boundaryMesh()[patchi]);

    celli_ = ppp.faceCells()[facei_];

    if (!ppp.parallel())
    {
        if (ppp.forwardT().size() == 1)
        {
            const tensor& T = ppp.forwardT()[0];
            transformPosition(T);
            static_cast<ParticleType&>(*this).transformProperties(T);
        }
        else
        {
            const tensor& T = ppp.forwardT()[facei_];
            transformPosition(T);
            static_cast<ParticleType&>(*this).transformProperties(T);
        }
    }
    else if (ppp.separated())
    {
        if (ppp.separation().size() == 1)
        {
            position_ -= ppp.separation()[0];
            static_cast<ParticleType&>(*this).transformProperties
            (
                -ppp.separation()[0]
            );
        }
        else
        {
            position_ -= ppp.separation()[facei_];
            static_cast<ParticleType&>(*this).transformProperties
            (
                -ppp.separation()[facei_]
            );
        }
    }

    // Reset the face index for the next tracking operation
    if (stepFraction_ > (1.0 - SMALL))
    {
        stepFraction_ = 1.0;
        facei_ = -1;
    }
    else
    {
        facei_ += ppp.start();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::Particle<ParticleType>::Particle
(
    const Cloud<ParticleType>& cloud,
    const vector& position,
    const label celli
)
:
    cloud_(cloud),
    position_(position),
    celli_(celli),
    facei_(-1),
    stepFraction_(0.0),
    origProc_(Pstream::myProcNo()),
    origId_(cloud_.getNewParticleID())
{}


template<class ParticleType>
Foam::Particle<ParticleType>::Particle(const Particle<ParticleType>& p)
:
    cloud_(p.cloud_),
    position_(p.position_),
    celli_(p.celli_),
    facei_(p.facei_),
    stepFraction_(p.stepFraction_),
    origProc_(p.origProc_),
    origId_(p.origId_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
template<class TrackData>
Foam::label Foam::Particle<ParticleType>::track
(
    const vector& endPosition,
    TrackData& td
)
{
    facei_ = -1;

    // Tracks to endPosition or stop on boundary
    while(!onBoundary() && stepFraction_ < 1.0 - SMALL)
    {
        stepFraction_ += trackToFace(endPosition, td)*(1.0 - stepFraction_);
    }

    return facei_;
}


template<class ParticleType>
Foam::label Foam::Particle<ParticleType>::track(const vector& endPosition)
{
    int dummyTd;
    return track(endPosition, dummyTd);
}


template<class ParticleType>
template<class TrackData>
Foam::scalar Foam::Particle<ParticleType>::trackToFace
(
    const vector& endPosition,
    TrackData& td
)
{
    scalar trackFraction = 0.0;

    if (cloud_.concaveCheck_)
    {
        if (cloud_.concaveCell()[celli_])
        {
            // Use a more careful tracking algorithm if the cell is concave
            trackToFaceConcave(trackFraction, endPosition, td);
        }
        else
        {
            // Use a more careful tracking algorithm if the cell is concave
            trackToFaceConvex(trackFraction, endPosition, td);
        }
    }
    else
    {
        // Use the original tracking algorithm if the cell is convex
        trackToFaceConvex(trackFraction, endPosition, td);
    }

    return trackFraction;
}


template<class ParticleType>
Foam::scalar Foam::Particle<ParticleType>::trackToFace
(
    const vector& endPosition
)
{
    int dummyTd;
    return trackToFace(endPosition, dummyTd);
}


template<class ParticleType>
void Foam::Particle<ParticleType>::transformPosition(const tensor& T)
{
    position_ = transform(T, position_);
}


template<class ParticleType>
void Foam::Particle<ParticleType>::transformProperties(const tensor&)
{}


template<class ParticleType>
void Foam::Particle<ParticleType>::transformProperties(const vector&)
{}


template<class ParticleType>
template<class TrackData>
bool Foam::Particle<ParticleType>::hitPatch
(
    const polyPatch&,
    TrackData&,
    const label
)
{
    return false;
}


template<class ParticleType>
template<class TrackData>
void Foam::Particle<ParticleType>::hitWedgePatch
(
    const wedgePolyPatch& wpp,
    TrackData&
)
{
    vector nf = wpp.faceAreas()[wpp.whichFace(facei_)];
    nf /= mag(nf);

    static_cast<ParticleType&>(*this).transformProperties(I - 2.0*nf*nf);
}


template<class ParticleType>
template<class TrackData>
void Foam::Particle<ParticleType>::hitSymmetryPatch
(
    const symmetryPolyPatch& spp,
    TrackData&
)
{
    vector nf = spp.faceAreas()[spp.whichFace(facei_)];
    nf /= mag(nf);

    static_cast<ParticleType&>(*this).transformProperties(I - 2.0*nf*nf);
}


template<class ParticleType>
template<class TrackData>
void Foam::Particle<ParticleType>::hitCyclicPatch
(
    const cyclicPolyPatch& cpp,
    TrackData&
)
{
    label patchFacei_ = cpp.whichFace(facei_);

    facei_ = cpp.transformGlobalFace(facei_);

    celli_ = cloud_.polyMesh_.faceOwner()[facei_];

    if (!cpp.parallel())
    {
        const tensor& T = cpp.transformT(patchFacei_);

        transformPosition(T);
        static_cast<ParticleType&>(*this).transformProperties(T);
    }
    else if (cpp.separated())
    {
        position_ += cpp.separation(patchFacei_);
        static_cast<ParticleType&>(*this).transformProperties
        (
            cpp.separation(patchFacei_)
        );
    }
}


template<class ParticleType>
template<class TrackData>
void Foam::Particle<ParticleType>::hitProcessorPatch
(
    const processorPolyPatch& spp,
    TrackData& td
)
{}


template<class ParticleType>
template<class TrackData>
void Foam::Particle<ParticleType>::hitWallPatch
(
    const wallPolyPatch& spp,
    TrackData&
)
{}


template<class ParticleType>
template<class TrackData>
void Foam::Particle<ParticleType>::hitPatch
(
    const polyPatch& spp,
    TrackData&
)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ParticleIO.C"

// ************************************************************************* //
