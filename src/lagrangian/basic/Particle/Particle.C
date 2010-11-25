/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "Particle.H"
#include "Cloud.H"
#include "wedgePolyPatch.H"
#include "symmetryPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "processorPolyPatch.H"
#include "transform.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ParticleType>
template<class TrackData>
void Foam::Particle<ParticleType>::prepareForParallelTransfer
(
    const label patchI,
    TrackData& td
)
{
    // Convert the face index to be local to the processor patch
    faceI_ = patchFace(patchI, faceI_);
}


template<class ParticleType>
template<class TrackData>
void Foam::Particle<ParticleType>::correctAfterParallelTransfer
(
    const label patchI,
    TrackData& td
)
{
    const processorPolyPatch& ppp =
        refCast<const processorPolyPatch>
        (cloud_.pMesh().boundaryMesh()[patchI]);

    cellI_ = ppp.faceCells()[faceI_];

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
            const tensor& T = ppp.forwardT()[faceI_];
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
            position_ -= ppp.separation()[faceI_];
            static_cast<ParticleType&>(*this).transformProperties
            (
                -ppp.separation()[faceI_]
            );
        }
    }

    tetFaceI_ = faceI_ + ppp.start();

    // Faces either side of a coupled patch have matched base indices,
    // tetPtI is specified relative to the base point, already and
    // opposite circulation directions by design, so if the vertices
    // are:
    // source:
    // face    (a b c d e f)
    // fPtI     0 1 2 3 4 5
    //            +
    // destination:
    // face    (a f e d c b)
    // fPtI     0 1 2 3 4 5
    //                  +
    // where a is the base point of the face are matching , and we
    // have fPtI = 1 on the source processor face, i.e. vertex b, then
    // this because of the face circulation direction change, vertex c
    // is the characterising point on the destination processor face,
    // giving the destination fPtI as:
    //     fPtI_d = f.size() - 1 - fPtI_s = 6 - 1 - 1 = 4
    // This relationship can be verified for other points and sizes of
    // face.

    tetPtI_ = cloud_.polyMesh_.faces()[tetFaceI_].size() - 1 - tetPtI_;

    // Reset the face index for the next tracking operation
    if (stepFraction_ > (1.0 - SMALL))
    {
        stepFraction_ = 1.0;
        faceI_ = -1;
    }
    else
    {
        faceI_ += ppp.start();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::Particle<ParticleType>::Particle
(
    const Cloud<ParticleType>& cloud,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI
)
:
    cloud_(cloud),
    position_(position),
    cellI_(cellI),
    faceI_(-1),
    stepFraction_(0.0),
    tetFaceI_(tetFaceI),
    tetPtI_(tetPtI),
    origProc_(Pstream::myProcNo()),
    origId_(cloud_.getNewParticleID())
{}


template<class ParticleType>
Foam::Particle<ParticleType>::Particle
(
    const Cloud<ParticleType>& cloud,
    const vector& position,
    const label cellI,
    bool doCellFacePt
)
:
    cloud_(cloud),
    position_(position),
    cellI_(cellI),
    faceI_(-1),
    stepFraction_(0.0),
    tetFaceI_(-1),
    tetPtI_(-1),
    origProc_(Pstream::myProcNo()),
    origId_(cloud_.getNewParticleID())
{
    if (doCellFacePt)
    {
        initCellFacePt();
    }
}


template<class ParticleType>
Foam::Particle<ParticleType>::Particle(const Particle<ParticleType>& p)
:
    cloud_(p.cloud_),
    position_(p.position_),
    cellI_(p.cellI_),
    faceI_(p.faceI_),
    stepFraction_(p.stepFraction_),
    tetFaceI_(p.tetFaceI_),
    tetPtI_(p.tetPtI_),
    origProc_(p.origProc_),
    origId_(p.origId_)
{}


template<class ParticleType>
Foam::Particle<ParticleType>::Particle
(
    const Particle<ParticleType>& p,
    const Cloud<ParticleType>& c
)
:
    cloud_(c),
    position_(p.position_),
    cellI_(p.cellI_),
    faceI_(p.faceI_),
    stepFraction_(p.stepFraction_),
    tetFaceI_(p.tetFaceI_),
    tetPtI_(p.tetPtI_),
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
    faceI_ = -1;

    // Tracks to endPosition or stop on boundary
    while (!onBoundary() && stepFraction_ < 1.0 - SMALL)
    {
        stepFraction_ += trackToFace(endPosition, td)*(1.0 - stepFraction_);
    }

    return faceI_;
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
    const polyMesh& mesh = cloud_.polyMesh_;

    const faceList& pFaces = mesh.faces();
    const pointField& pPts = mesh.points();
    const vectorField& pC = mesh.cellCentres();

    faceI_ = -1;

    // Pout<< "Particle " << origId_ << " " << origProc_
    //     << " Tracking from " << position_
    //     << " to " << endPosition
    //     << endl;

    // Pout<< "stepFraction " << stepFraction_ << nl
    //     << "cellI " << cellI_ << nl
    //     << "tetFaceI " << tetFaceI_ << nl
    //     << "tetPtI " << tetPtI_
    //     << endl;

    scalar trackFraction = 0.0;

    // Minimum tetrahedron decomposition of each cell of the mesh into
    // using the cell centre, base point on face, and further two
    // points on the face.  For each face of n points, there are n - 2
    // tets generated.

    // The points for each tet are organised to match those used in the
    // tetrahedron class, supplying them in the order:
    //     Cc, basePt, pA, pB
    // where:
    //   + Cc is the cell centre;
    //   + basePt is the base point on the face;
    //   + pA and pB are the remaining points on the face, such that
    //     the circulation, {basePt, pA, pB} produces a positive
    //     normal by the right-hand rule.  pA and pB are chosen from
    //     tetPtI_ do accomplish this depending if the cell owns the
    //     face, tetPtI_ is the vertex that characterises the tet, and
    //     is the first vertex on the tet when circulating around the
    //     face. Therefore, the same tetPtI represents the same face
    //     triangle for both the owner and neighbour cell.
    //
    // Each tet has its four triangles represented in the same order:
    // 0) tri joining a tet to the tet across the face in next cell.
    //    This is the triangle opposite Cc.
    // 1) tri joining a tet to the tet that is in the same cell, but
    //    belongs to the face that shares the edge of the current face
    //    that doesn't contain basePt.  This is the triangle opposite
    //    basePt.

    // 2) tri joining a tet to the tet that is in the same cell, but
    //    belongs to the face that shares the tet-edge (basePt - pB).
    //    This may be on the same face, or a different one.  This is
    //    the triangle opposite basePt.  This is the triangle opposite
    //    pA.

    // 4) tri joining a tet to the tet that is in the same cell, but
    //    belongs to the face that shares the tet-edge (basePt - pA).
    //    This may be on the same face, or a different one.  This is
    //    the triangle opposite basePt.  This is the triangle opposite
    //    pA.

    // Which tri (0..3) of the tet has been crossed
    label triI = -1;

    // Determine which face was actually crossed.  lambdaMin < SMALL
    // is considered a trigger for a tracking correction towards the
    // current tet centre.
    scalar lambdaMin = VGREAT;

    DynamicList<label>& tris = cloud_.labels_;

    // Tet indices that will be set by hitWallFaces if a wall face is
    // to be hit, or are set when any wall tri of a tet is hit.
    // Carries the description of the tet on which the cell face has
    // been hit.  For the case of being set in hitWallFaces, this may
    // be a different tet to the one that the particle occupies.
    tetIndices faceHitTetIs;

    do
    {
        if (triI != -1)
        {
            // Change tet ownership because a tri face has been crossed
            tetNeighbour(triI);
        }

        const Foam::face& f = pFaces[tetFaceI_];

        bool own = (mesh.faceOwner()[tetFaceI_] == cellI_);

        label tetBasePtI = mesh.tetBasePtIs()[tetFaceI_];

        label basePtI = f[tetBasePtI];

        label facePtI = (tetPtI_ + tetBasePtI) % f.size();
        label otherFacePtI = f.fcIndex(facePtI);

        label fPtAI = -1;
        label fPtBI = -1;

        if (own)
        {
            fPtAI = facePtI;
            fPtBI = otherFacePtI;
        }
        else
        {
            fPtAI = otherFacePtI;
            fPtBI = facePtI;
        }

        tetPointRef tet
        (
            pC[cellI_],
            pPts[basePtI],
            pPts[f[fPtAI]],
            pPts[f[fPtBI]]
        );

        if (lambdaMin < SMALL)
        {
            // Apply tracking correction towards tet centre

            position_ +=
                Cloud<ParticleType>::trackingCorrectionTol
               *(tet.centre() - position_);

            cloud_.trackingRescue();

            return trackFraction;
        }

        if (triI != -1 && mesh.moving())
        {
            // Mesh motion requires stepFraction to be correct for
            // each tracking portion, so trackToFace must return after
            // every lambda calculation.
            return trackFraction;
        }

        FixedList<vector, 4> tetAreas;

        tetAreas[0] = tet.Sa();
        tetAreas[1] = tet.Sb();
        tetAreas[2] = tet.Sc();
        tetAreas[3] = tet.Sd();

        FixedList<label, 4> tetPlaneBasePtIs;

        tetPlaneBasePtIs[0] = basePtI;
        tetPlaneBasePtIs[1] = f[fPtAI];
        tetPlaneBasePtIs[2] = basePtI;
        tetPlaneBasePtIs[3] = basePtI;

        findTris(endPosition, tris, tet, tetAreas, tetPlaneBasePtIs);

        // Reset variables for new track
        triI = -1;
        lambdaMin = VGREAT;

        // Pout<< "tris " << tris << endl;

        // Sets a value for lambdaMin and faceI_ if a wall face is hit
        // by the track.
        hitWallFaces(position_, endPosition, lambdaMin, faceHitTetIs);

        // Did not hit any tet tri faces, and no wall face has been
        // found to hit.
        if (tris.empty() && faceI_ < 0)
        {
            position_ = endPosition;

            return 1.0;
        }
        else
        {
            // Loop over all found tris and see if any of them find a
            // lambda value smaller than that found for a wall face.
            forAll(tris, i)
            {
                label tI = tris[i];

                scalar lam = tetLambda
                (
                    position_,
                    endPosition,
                    triI,
                    tetAreas[tI],
                    tetPlaneBasePtIs[tI],
                    cellI_,
                    tetFaceI_,
                    tetPtI_
                );

                if (lam < lambdaMin)
                {
                    lambdaMin = lam;

                    triI = tI;
                }
            }
        }

        if (triI == 0)
        {
            // This must be a cell face crossing
            faceI_ = tetFaceI_;

            // Set the faceHitTetIs to those for the current tet in case a
            // wall interaction is required with the cell face
            faceHitTetIs = tetIndices
            (
                cellI_,
                tetFaceI_,
                tetBasePtI,
                fPtAI,
                fPtBI,
                tetPtI_
            );
        }
        else if (triI > 0)
        {
            // A tri was found to be crossed before a wall face was hit (if any)
            faceI_ = -1;
        }

        // Pout<< "track loop " << position_ << " " << endPosition << nl
        //     << "    " << cellI_
        //     << "    " << faceI_
        //     << " " << tetFaceI_
        //     << " " << tetPtI_
        //     << " " << triI
        //     << " " << lambdaMin
        //     << " " << trackFraction
        //     << endl;

        // Pout<< "# Tracking loop tet "
        //     << origId_ << " " << origProc_<< nl
        //     << "# face: " << tetFaceI_ << nl
        //     << "# tetPtI: " << tetPtI_ << nl
        //     << "# tetBasePtI: " << mesh.tetBasePtIs()[tetFaceI_] << nl
        //     << "# tet.mag(): " << tet.mag() << nl
        //     << "# tet.quality(): " << tet.quality()
        //     << endl;

        // meshTools::writeOBJ(Pout, tet.a());
        // meshTools::writeOBJ(Pout, tet.b());
        // meshTools::writeOBJ(Pout, tet.c());
        // meshTools::writeOBJ(Pout, tet.d());

        // Pout<< "f 1 3 2" << nl
        //     << "f 2 3 4" << nl
        //     << "f 1 4 3" << nl
        //     << "f 1 2 4" << endl;

        // The particle can be 'outside' the tet.  This will yield a
        // lambda larger than 1, or smaller than 0.  For values < 0,
        // the particle travels away from the tet and we don't move
        // the particle, only change tet/cell.  For values larger than
        // 1, we move the particle to endPosition before the tet/cell
        // change.
        if (lambdaMin > SMALL)
        {
            if (lambdaMin <= 1.0)
            {
                trackFraction += lambdaMin*(1 - trackFraction);

                position_ += lambdaMin*(endPosition - position_);
            }
            else
            {
                trackFraction = 1.0;

                position_ = endPosition;
            }
        }
        else
        {
            // Set lambdaMin to zero to force a towards-tet-centre
            // correction.
            lambdaMin = 0.0;
        }

    } while (faceI_ < 0);

    ParticleType& p = static_cast<ParticleType&>(*this);
    p.hitFace(td);

    if (cloud_.internalFace(faceI_))
    {
        // Change tet ownership because a tri face has been crossed,
        // in general this is:
        //     tetNeighbour(triI);
        // but triI must be 0;
        // No modifications are required for triI = 0, no call required to
        //     tetNeighbour(0);

        if (cellI_ == mesh.faceOwner()[faceI_])
        {
            cellI_ = mesh.faceNeighbour()[faceI_];
        }
        else if (cellI_ == mesh.faceNeighbour()[faceI_])
        {
            cellI_ = mesh.faceOwner()[faceI_];
        }
        else
        {
            FatalErrorIn
            (
                "Particle::trackToFace(const vector&, TrackData&)"
            )   << "addressing failure" << nl
                << abort(FatalError);
        }
    }
    else
    {
        label origFaceI = faceI_;
        label patchI = patch(faceI_);

        // No action taken for tetPtI_ for tetFaceI_ here, handled by
        // patch interaction call or later during processor transfer.

        if
        (
            !p.hitPatch
            (
                mesh.boundaryMesh()[patchI],
                td,
                patchI,
                trackFraction,
                faceHitTetIs
            )
        )
        {
            // Did patch interaction model switch patches?
            if (faceI_ != origFaceI)
            {
                patchI = patch(faceI_);
            }

            const polyPatch& patch = mesh.boundaryMesh()[patchI];

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
                    static_cast<const wallPolyPatch&>(patch), td, faceHitTetIs
                );
            }
            else
            {
                p.hitPatch(patch, td);
            }
        }
    }

    if (lambdaMin < SMALL)
    {
        // Apply tracking correction towards tet centre.
        // Generate current tet to find centre to apply correction.

        tetPointRef tet = currentTet();

        position_ +=
            Cloud<ParticleType>::trackingCorrectionTol
           *(tet.centre() - position_);

        if
        (
            cloud_.hasWallImpactDistance()
            && !cloud_.internalFace(faceHitTetIs.face())
            && cloud_.cellHasWallFaces()[faceHitTetIs.cell()]
        )
        {
            const polyBoundaryMesh& patches = mesh.boundaryMesh();

            label fI = faceHitTetIs.face();

            label patchI = patches.patchID()[fI - mesh.nInternalFaces()];

            if (isA<wallPolyPatch>(patches[patchI]))
            {
                // In the case of collision with a wall where there is
                // a non-zero wallImpactDistance, it is possible for
                // there to be a tracking correction required to bring
                // the particle into the domain, but the position of
                // the particle is further from the wall than the tet
                // centre, in which case the normal correction can be
                // counter-productive, i.e. pushes the particle
                // further out of the domain.  In this case it is the
                // position that hit the wall that is in need of a
                // rescue correction.

                triPointRef wallTri = faceHitTetIs.faceTri(mesh);

                tetPointRef wallTet = faceHitTetIs.tet(mesh);

                vector nHat = wallTri.normal();
                nHat /= mag(nHat);

                const scalar r = p.wallImpactDistance(nHat);

                // Removing (approximately) the wallTri normal
                // component of the existing correction, to avoid the
                // situation where the existing correction in the wall
                // normal direction is larger towards the wall than
                // the new correction is away from it.
                position_ +=
                    Cloud<ParticleType>::trackingCorrectionTol
                   *(
                        (wallTet.centre() - (position_ + r*nHat))
                      - (nHat & (tet.centre() - position_))*nHat
                    );
            }
        }

        cloud_.trackingRescue();
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
void Foam::Particle<ParticleType>::hitFace(TrackData&)
{}


template<class ParticleType>
template<class TrackData>
bool Foam::Particle<ParticleType>::hitPatch
(
    const polyPatch&,
    TrackData&,
    const label,
    const scalar,
    const tetIndices&
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
    FatalErrorIn
    (
        "void Foam::Particle<ParticleType>::hitWedgePatch"
        "("
            "const wedgePolyPatch& wpp, "
            "TrackData&"
        ")"
    )   << "Hitting a wedge patch should not be possible."
        << abort(FatalError);

    vector nf = normal();
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
    vector nf = normal();
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
    // label patchFaceI_ = cpp.whichFace(faceI_);

    faceI_ = cpp.transformGlobalFace(faceI_);

    cellI_ = cloud_.polyMesh_.faceOwner()[faceI_];

    tetFaceI_ = faceI_;

    // See note in correctAfterParallelTransfer for tetPtI_ addressing.
    tetPtI_ = cloud_.polyMesh_.faces()[tetFaceI_].size() - 1 - tetPtI_;

    // Now the particle is on the receiving side

    if (!cpp.parallel())
    {
        const tensor& T = cpp.reverseT()[0];

        transformPosition(T);
        static_cast<ParticleType&>(*this).transformProperties(T);
    }
    else if (cpp.separated())
    {
        position_ += cpp.separation()[0];
        static_cast<ParticleType&>(*this).transformProperties
        (
            cpp.separation()[0]
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
    TrackData&,
    const tetIndices&
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


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

template<class ParticleType>
bool Foam::operator==
(
    const Particle<ParticleType>& pA,
    const Particle<ParticleType>& pB
)
{
    return
    (
        pA.origProc() == pB.origProc()
     && pA.origId() == pB.origId()
    );
}


template<class ParticleType>
bool Foam::operator!=
(
    const Particle<ParticleType>& pA,
    const Particle<ParticleType>& pB
)
{
    return !(pA == pB);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ParticleIO.C"

// ************************************************************************* //
