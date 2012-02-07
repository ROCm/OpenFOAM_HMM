/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

\*----------------------------------------------------------------------------*/

#include "wallBoundedStreamLineParticle.H"
#include "vectorFieldIOField.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
//    defineParticleTypeNameAndDebug(wallBoundedStreamLineParticle, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//// Check position is inside tet
//void Foam::wallBoundedStreamLineParticle::checkInside() const
//{
//    const tetIndices ti(currentTetIndices());
//    const tetPointRef tpr(ti.tet(mesh_));
//    if (!tpr.inside(position()))
//    {
//        FatalErrorIn("wallBoundedStreamLineParticle::checkInside(..)")
//            << "Particle:" //<< static_cast<const particle&>(*this)
//            << info()
//            << "is not inside " << tpr
//            << abort(FatalError);
//    }
//}
//
//
//void Foam::wallBoundedStreamLineParticle::checkOnEdge() const
//{
//    // Check that edge (as indicated by meshEdgeStart_, diagEdge_) is
//    // indeed one that contains the position.
//    const edge e = currentEdge();
//
//    linePointRef ln(e.line(mesh_.points()));
//
//    pointHit ph(ln.nearestDist(position()));
//
//    if (ph.distance() > 1E-6)
//    {
//        FatalErrorIn
//        (
//            "wallBoundedStreamLineParticle::checkOnEdge()"
//        )   << "Problem :"
//            << " particle:" //<< static_cast<const particle&>(*this)
//            << info()
//            << "edge:" << e
//            << " at:" << ln
//            << " distance:" << ph.distance()
//            << abort(FatalError);
//    }
//}
//
//
//void Foam::wallBoundedStreamLineParticle::checkOnTriangle(const point& p)
//const
//{
//    const triFace tri(currentTetIndices().faceTriIs(mesh_));
//    pointHit ph = tri.nearestPoint(p, mesh_.points());
//    if (ph.distance() > 1e-9)
//    {
//        FatalErrorIn
//        (
//            "wallBoundedStreamLineParticle::checkOnTriangle(const point&)"
//        )   << "Problem :"
//            << " particle:" //<< static_cast<const particle&>(*this)
//            << info()
//            << "point:" << p
//            << " distance:" << ph.distance()
//            << abort(FatalError);
//    }
//}


// Construct the edge the particle is on (according to meshEdgeStart_,
// diagEdge_)
Foam::edge Foam::wallBoundedStreamLineParticle::currentEdge() const
{
    if ((meshEdgeStart_ != -1) == (diagEdge_ != -1))
    {
        FatalErrorIn("wallBoundedStreamLineParticle::currentEdge() const")
            << "Particle:" //<< static_cast<const particle&>(*this)
            << info()
            << "cannot both be on a mesh edge and a face-diagonal edge."
            << " meshEdgeStart_:" << meshEdgeStart_
            << " diagEdge_:" << diagEdge_
            << abort(FatalError);
    }

    const Foam::face& f = mesh_.faces()[tetFace()];

    if (meshEdgeStart_ != -1)
    {
        return edge(f[meshEdgeStart_], f.nextLabel(meshEdgeStart_));
    }
    else
    {
        label faceBasePtI = mesh_.tetBasePtIs()[tetFace()];
        label diagPtI = (faceBasePtI+diagEdge_)%f.size();
        return edge(f[faceBasePtI], f[diagPtI]);
    }
}


void Foam::wallBoundedStreamLineParticle::crossEdgeConnectedFace
(
    const edge& meshEdge
)
{
    //label oldFaceI = tetFace();

    // Update tetFace, tetPt
    particle::crossEdgeConnectedFace(cell(), tetFace(), tetPt(), meshEdge);

    // Update face to be same as tracking one
    face() = tetFace();


    // And adapt meshEdgeStart_.
    const Foam::face& f = mesh_.faces()[tetFace()];
    label fp = findIndex(f, meshEdge[0]);

    if (f.nextLabel(fp) == meshEdge[1])
    {
        meshEdgeStart_ = fp;
    }
    else
    {
        label fpMin1 = f.rcIndex(fp);

        if (f[fpMin1] == meshEdge[1])
        {
            meshEdgeStart_ = fpMin1;
        }
        else
        {
            FatalErrorIn
            (
                "wallBoundedStreamLineParticle::crossEdgeConnectedFace"
                "(const edge&)"
            )   << "Problem :"
                << " particle:" //<< static_cast<const particle&>(*this)
                << info()
                << "face:" << tetFace()
                << " verts:" << f
                << " meshEdge:" << meshEdge
                << abort(FatalError);
        }
    }

    diagEdge_ = -1;

    //Pout<< "    crossed meshEdge "
    //    << meshEdge.line(mesh().points())
    //    << " from face:" << oldFaceI
    //    << " to face:" << tetFace() << endl;


    // Check that still on same mesh edge

    const edge eNew(f[meshEdgeStart_], f.nextLabel(meshEdgeStart_));
    if (eNew != meshEdge)
    {
        FatalErrorIn
        (
            "wallBoundedStreamLineParticle::crossEdgeConnectedFace"
            "(const edge&)"
        )   << "Problem" << abort(FatalError);
    }


    // Check that edge (as indicated by meshEdgeStart_) is indeed one that
    // contains the position.
    //checkOnEdge();
}


void Foam::wallBoundedStreamLineParticle::crossDiagonalEdge()
{
    if (diagEdge_ == -1)
    {
        FatalErrorIn("wallBoundedStreamLineParticle::crossDiagonalEdge()")
            << "Particle:" //<< static_cast<const particle&>(*this)
            << info()
            << "not on a diagonal edge" << abort(FatalError);
    }
    if (meshEdgeStart_ != -1)
    {
        FatalErrorIn("wallBoundedStreamLineParticle::crossDiagonalEdge()")
            << "Particle:" //<< static_cast<const particle&>(*this)
            << info()
            << "meshEdgeStart_:" << meshEdgeStart_ << abort(FatalError);
    }

    //label oldTetPt = tetPt();

    const Foam::face& f = mesh_.faces()[tetFace()];

    // tetPtI starts from 1, goes up to f.size()-2

    if (tetPt() == diagEdge_)
    {
        tetPt() = f.rcIndex(tetPt());
    }
    else
    {
        label nextTetPt = f.fcIndex(tetPt());
        if (diagEdge_ == nextTetPt)
        {
            tetPt() = nextTetPt;
        }
        else
        {
            FatalErrorIn("wallBoundedStreamLineParticle::crossDiagonalEdge()")
                << "Particle:" //<< static_cast<const particle&>(*this)
                << info()
                << "tetPt:" << tetPt()
                << " diagEdge:" << diagEdge_ << abort(FatalError);
        }
    }

    meshEdgeStart_ = -1;

    //Pout<< "    crossed diagonalEdge "
    //    << currentEdge().line(mesh().points())
    //    << " from tetPt:" << oldTetPt
    //    << " to tetPt:" << tetPt() << endl;
}


//- Track through a single triangle.
// Gets passed tet+triangle the particle is in. Updates position() but nothing
// else. Returns the triangle edge the particle is now on.
Foam::scalar Foam::wallBoundedStreamLineParticle::trackFaceTri
(
    const vector& endPosition,
    label& minEdgeI
)
{
    // Track p from position to endPosition
    const triFace tri(currentTetIndices().faceTriIs(mesh_));
    vector n = tri.normal(mesh_.points());
    //if (mag(n) < sqr(SMALL))
    //{
    //    FatalErrorIn("wallBoundedStreamLineParticle::trackFaceTri(..)")
    //        << "Small triangle." //<< static_cast<const particle&>(*this)
    //        << info()
    //        << "n:" << n
    //        << abort(FatalError);
    //}
    n /= mag(n)+VSMALL;

    // Check which edge intersects the trajectory.
    // Project trajectory onto triangle
    minEdgeI = -1;
    scalar minS = 1;        // end position

    //const point oldPosition(position());


    edge currentE(-1, -1);
    if (meshEdgeStart_ != -1 || diagEdge_ != -1)
    {
        currentE = currentEdge();
    }

    // Determine path along line position+s*d to see where intersections
    // are.

    forAll(tri, i)
    {
        label j = tri.fcIndex(i);

        const point& pt0 = mesh_.points()[tri[i]];
        const point& pt1 = mesh_.points()[tri[j]];

        if (edge(tri[i], tri[j]) == currentE)
        {
            // Do not check particle is on
            continue;
        }

        // Outwards pointing normal
        vector edgeNormal = (pt1-pt0)^n;

        //if (mag(edgeNormal) < SMALL)
        //{
        //    FatalErrorIn("wallBoundedStreamLineParticle::trackFaceTri(..)")
        //        << "Edge not perpendicular to triangle."
        //        //<< static_cast<const particle&>(*this)
        //        << info()
        //        << "triangle n:" << n
        //        << " edgeNormal:" << edgeNormal
        //        << " on tri:" << tri
        //        << " at:" << pt0
        //        << " at:" << pt1
        //        << abort(FatalError);
        //}


        edgeNormal /= mag(edgeNormal)+VSMALL;

        // Determine whether position and end point on either side of edge.
        scalar sEnd = (endPosition-pt0)&edgeNormal;
        if (sEnd >= 0)
        {
            // endPos is outside triangle. position() should always be
            // inside.
            scalar sStart = (position()-pt0)&edgeNormal;
            if (mag(sEnd-sStart) > VSMALL)
            {
                scalar s = sStart/(sStart-sEnd);

                if (s >= 0 && s < minS)
                {
                    minS = s;
                    minEdgeI = i;
                }
            }
        }
    }

    if (minEdgeI != -1)
    {
        position() += minS*(endPosition-position());
    }
    else
    {
        // Did not hit any edge so tracked to the end position. Set position
        // without any calculation to avoid truncation errors.
        position() = endPosition;
        minS = 1.0;
    }

    // Project position() back onto plane of triangle
    const point& triPt = mesh_.points()[tri[0]];
    position() -= ((position()-triPt)&n)*n;


    //Pout<< "    tracked from:" << oldPosition << " to:" << position()
    //    << " projectedEnd:" << endPosition
    //    << " at s:" << minS << endl;
    //if (minEdgeI != -1)
    //{
    //    Pout<< "    on edge:" << minEdgeI
    //        << " on edge:"
    //        << mesh_.points()[tri[minEdgeI]]
    //        << mesh_.points()[tri[tri.fcIndex(minEdgeI)]]
    //        << endl;
    //}

    return minS;
}


// See if the current triangle has got a point on the
// correct side of the edge.
bool Foam::wallBoundedStreamLineParticle::isTriAlongTrack
(
    const point& endPosition
) const
{
    const triFace triVerts(currentTetIndices().faceTriIs(mesh_));
    const edge currentE = currentEdge();

    //if (debug)
    {
        if
        (
            currentE[0] == currentE[1]
         || findIndex(triVerts, currentE[0]) == -1
         || findIndex(triVerts, currentE[1]) == -1
        )
        {
            FatalErrorIn
            (
                "wallBoundedStreamLineParticle::isTriAlongTrack"
                "(const point&)"
            )   << "Edge " << currentE << " not on triangle " << triVerts
                << info()
                << abort(FatalError);
        }
    }


    const vector dir = endPosition-position();

    // Get normal of currentE
    vector n = triVerts.normal(mesh_.points());
    n /= mag(n);

    forAll(triVerts, i)
    {
        label j = triVerts.fcIndex(i);
        const point& pt0 = mesh_.points()[triVerts[i]];
        const point& pt1 = mesh_.points()[triVerts[j]];

        if (edge(triVerts[i], triVerts[j]) == currentE)
        {
            vector edgeNormal = (pt1-pt0)^n;
            return (dir&edgeNormal) < 0;
        }
    }

    FatalErrorIn
    (
        "wallBoundedStreamLineParticle::isTriAlongTrack"
        "(const point&)"
    )   << "Problem" << abort(FatalError);

    return false;
}


void Foam::wallBoundedStreamLineParticle::patchInteraction
(
    trackingData& td,
    const scalar trackFraction
)
{
    typedef trackingData::cloudType cloudType;
    typedef cloudType::particleType particleType;

    particleType& p = static_cast<particleType&>(*this);
    p.hitFace(td);

    if (!internalFace(faceI_))
    {
        label origFaceI = faceI_;
        label patchI = patch(faceI_);

        // No action taken for tetPtI_ for tetFaceI_ here, handled by
        // patch interaction call or later during processor transfer.


        // Dummy tet indices. What to do here?
        tetIndices faceHitTetIs;

        if
        (
            !p.hitPatch
            (
                mesh_.boundaryMesh()[patchI],
                td,
                patchI,
                trackFraction,
                faceHitTetIs
            )
        )
        {
            // Did patch interaction model switch patches?
            // Note: recalculate meshEdgeStart_, diagEdge_!
            if (faceI_ != origFaceI)
            {
                patchI = patch(faceI_);
            }

            const polyPatch& patch = mesh_.boundaryMesh()[patchI];

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
}

//- Track particle to a given position and returns 1.0 if the
//  trajectory is completed without hitting a face otherwise
//  stops at the face and returns the fraction of the trajectory
//  completed.
//  on entry 'stepFraction()' should be set to the fraction of the
//  time-step at which the tracking starts.
Foam::scalar Foam::wallBoundedStreamLineParticle::trackToEdge
(
    trackingData& td,
    const vector& endPosition
)
{
    // Are we on a track face? If not we do a topological walk.

    // Particle:
    // - cell_              always set
    // - tetFace_, tetPt_   always set (these identify tet particle is in)
    // - optionally meshEdgeStart_ or  diagEdge_ set (edge particle is on)

    //checkInside();
    //checkOnTriangle(position());
    //if (meshEdgeStart_ != -1 || diagEdge_ != -1)
    //{
    //    checkOnEdge();
    //}

    scalar trackFraction = 0.0;

    if (!td.isWallPatch_[tetFace()])
    {
        // Don't track across face. Just walk in cell. Particle is on
        // mesh edge (as indicated by meshEdgeStart_).
        const edge meshEdge(currentEdge());

        // If internal face check whether to go to neighbour cell or just
        // check to the other internal tet on the edge.
        if (mesh_.isInternalFace(tetFace()))
        {
            label nbrCellI =
            (
                cellI_ == mesh_.faceOwner()[faceI_]
              ? mesh_.faceNeighbour()[faceI_]
              : mesh_.faceOwner()[faceI_]
            );
            // Check angle to nbrCell tet. Is it in the direction of the
            // endposition? I.e. since volume of nbr tet is positive the
            // tracking direction should be into the tet.
            tetIndices nbrTi(nbrCellI, tetFaceI_, tetPtI_, mesh_);
            if ((nbrTi.faceTri(mesh_).normal() & (endPosition-position())) < 0)
            {
                // Change into nbrCell. No need to change tetFace, tetPt.
                //Pout<< "    crossed from cell:" << cellI_
                //    << " into " << nbrCellI << endl;
                cellI_ = nbrCellI;
                patchInteraction(td, trackFraction);
            }
            else
            {
                // Walk to other face on edge. Changes tetFace, tetPt but not
                // cell.
                crossEdgeConnectedFace(meshEdge);
                patchInteraction(td, trackFraction);
            }
        }
        else
        {
            // Walk to other face on edge. This might give loop since
            // particle should have been removed?
            crossEdgeConnectedFace(meshEdge);
            patchInteraction(td, trackFraction);
        }
    }
    else
    {
        // We're inside a tet on the wall. Check if the current tet is
        // the one to cross. If not we cross into the neighbouring triangle.

        if (mesh_.isInternalFace(tetFace()))
        {
            FatalErrorIn
            (
                "wallBoundedStreamLineParticle::trackToEdge"
                "(trackingData&, const vector&)"
            )   << "Can only track on boundary faces."
                << " Face:" << tetFace()
                << " at:" << mesh_.faceCentres()[tetFace()]
                << abort(FatalError);
        }

        point projectedEndPosition = endPosition;
        // Remove normal component
        {
            const triFace tri(currentTetIndices().faceTriIs(mesh_));
            vector n = tri.normal(mesh_.points());
            n /= mag(n);
            const point& basePt = mesh_.points()[tri[0]];
            projectedEndPosition -= ((projectedEndPosition-basePt)&n)*n;
        }


        bool doTrack = false;
        if (meshEdgeStart_ == -1 && diagEdge_ == -1)
        {
            // We're starting and not yet on an edge.
            doTrack = true;
        }
        else
        {
            // See if the current triangle has got a point on the
            // correct side of the edge.
            doTrack = isTriAlongTrack(projectedEndPosition);
        }


        if (doTrack)
        {
            // Track across triangle. Return triangle edge crossed.
            label triEdgeI = -1;
            trackFraction = trackFaceTri(projectedEndPosition, triEdgeI);

            if (triEdgeI == -1)
            {
                // Reached endpoint
                //checkInside();
                diagEdge_ = -1;
                meshEdgeStart_ = -1;
                return trackFraction;
            }

            const tetIndices ti(currentTetIndices());

            // Triangle (faceTriIs) gets constructed from
            //    f[faceBasePtI_],
            //    f[facePtAI_],
            //    f[facePtBI_]
            //
            // So edge indices are:
            // 0 : edge between faceBasePtI_ and facePtAI_
            // 1 : edge between facePtAI_ and facePtBI_ (is always a real edge)
            // 2 : edge between facePtBI_ and faceBasePtI_

            const Foam::face& f = mesh_.faces()[ti.face()];
            const label fp0 = ti.faceBasePt();

            if (triEdgeI == 0)
            {
                if (ti.facePtA() == f.fcIndex(fp0))
                {
                    //Pout<< "Real edge." << endl;
                    diagEdge_ = -1;
                    meshEdgeStart_ = fp0;
                    //checkOnEdge();
                    crossEdgeConnectedFace(currentEdge());
                    patchInteraction(td, trackFraction);
                }
                else if (ti.facePtA() == f.rcIndex(fp0))
                {
                    //Note: should not happen since boundary face so owner
                    //Pout<< "Real edge." << endl;
                    FatalErrorIn("shold not happend") << info()
                        << abort(FatalError);

                    diagEdge_ = -1;
                    meshEdgeStart_ = f.rcIndex(fp0);
                    //checkOnEdge();
                    crossEdgeConnectedFace(currentEdge());
                    patchInteraction(td, trackFraction);
                }
                else
                {
                    // Get index of triangle on other side of edge.
                    diagEdge_ = ti.facePtA()-fp0;
                    if (diagEdge_ < 0)
                    {
                        diagEdge_ += f.size();
                    }
                    meshEdgeStart_ = -1;
                    //checkOnEdge();
                    crossDiagonalEdge();
                }
            }
            else if (triEdgeI == 1)
            {
                //Pout<< "Real edge." << endl;
                diagEdge_ = -1;
                meshEdgeStart_ = ti.facePtA();
                //checkOnEdge();
                crossEdgeConnectedFace(currentEdge());
                patchInteraction(td, trackFraction);
            }
            else // if (triEdgeI == 2)
            {
                if (ti.facePtB() == f.rcIndex(fp0))
                {
                    //Pout<< "Real edge." << endl;
                    diagEdge_ = -1;
                    meshEdgeStart_ = ti.facePtB();
                    //checkOnEdge();
                    crossEdgeConnectedFace(currentEdge());
                    patchInteraction(td, trackFraction);
                }
                else if (ti.facePtB() == f.fcIndex(fp0))
                {
                    //Note: should not happen since boundary face so owner
                    //Pout<< "Real edge." << endl;
                    FatalErrorIn("shold not happend") << info()
                        << abort(FatalError);

                    diagEdge_ = -1;
                    meshEdgeStart_ = fp0;
                    //checkOnEdge();
                    crossEdgeConnectedFace(currentEdge());
                    patchInteraction(td, trackFraction);
                }
                else
                {
                    //Pout<< "Triangle edge." << endl;
                    // Get index of triangle on other side of edge.
                    diagEdge_ = ti.facePtB()-fp0;
                    if (diagEdge_ < 0)
                    {
                        diagEdge_ += f.size();
                    }
                    meshEdgeStart_ = -1;
                    //checkOnEdge();
                    crossDiagonalEdge();
                }
            }
        }
        else
        {
            // Current tet is not the right one. Check the neighbour tet.

            if (meshEdgeStart_ != -1)
            {
                // Particle is on mesh edge so change into other face on cell
                crossEdgeConnectedFace(currentEdge());
                //checkOnEdge();
                patchInteraction(td, trackFraction);
            }
            else
            {
                // Particle is on diagonal edge so change into the other
                // triangle.
                crossDiagonalEdge();
                //checkOnEdge();
            }
        }
    }

    //checkInside();

    return trackFraction;
}


Foam::vector Foam::wallBoundedStreamLineParticle::interpolateFields
(
    const trackingData& td,
    const point& position,
    const label cellI,
    const label faceI
)
{
    if (cellI == -1)
    {
        FatalErrorIn("wallBoundedStreamLineParticle::interpolateFields(..)")
            << "Cell:" << cellI << abort(FatalError);
    }

    const tetIndices ti = currentTetIndices();

    const vector U = td.vvInterp_[td.UIndex_].interpolate
    (
        position,
        ti,     //cellI,
        faceI
    );

    // Check if at different position
    if
    (
       !sampledPositions_.size()
     || magSqr(sampledPositions_.last()-position) > Foam::sqr(SMALL)
    )
    {
        // Store the location
        sampledPositions_.append(position);

        // Store the scalar fields
        sampledScalars_.setSize(td.vsInterp_.size());
        forAll(td.vsInterp_, scalarI)
        {
            sampledScalars_[scalarI].append
            (
                td.vsInterp_[scalarI].interpolate
                (
                    position,
                    ti,     //cellI,
                    faceI
                )
            );
        }

        // Store the vector fields
        sampledVectors_.setSize(td.vvInterp_.size());
        forAll(td.vvInterp_, vectorI)
        {
            vector positionU;
            if (vectorI == td.UIndex_)
            {
                positionU = U;
            }
            else
            {
                positionU = td.vvInterp_[vectorI].interpolate
                (
                    position,
                    ti,     //cellI,
                    faceI
                );
            }
            sampledVectors_[vectorI].append(positionU);
        }
    }

    return U;
}


Foam::vector Foam::wallBoundedStreamLineParticle::sample
(
    trackingData& td
)
{
    vector U = interpolateFields(td, position(), cell(), tetFace());

    if (!td.trackForward_)
    {
        U = -U;
    }

    scalar magU = mag(U);

    if (magU < SMALL)
    {
        // Stagnant particle. Might as well stop
        lifeTime_ = 0;
    }

    U /= magU;

    return U;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoundedStreamLineParticle::wallBoundedStreamLineParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const label meshEdgeStart,
    const label diagEdge,
    const label lifeTime
)
:
    particle(mesh, position, cellI, tetFaceI, tetPtI),
    meshEdgeStart_(meshEdgeStart),
    diagEdge_(diagEdge),
    lifeTime_(lifeTime)
{
    //checkInside();

    //if (meshEdgeStart_ != -1 || diagEdge_ != -1)
    //{
    //    checkOnEdge();
    //}

    // Unfortunately have no access to trackdata so cannot check if particle
    // is on a wallPatch or has an mesh edge set (either of which is
    // a requirement).
}


Foam::wallBoundedStreamLineParticle::wallBoundedStreamLineParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    particle(mesh, is, readFields)
{
    if (readFields)
    {
        List<scalarList> sampledScalars;
        List<vectorList> sampledVectors;

        is  >> meshEdgeStart_ >> diagEdge_ >> lifeTime_
            >> sampledPositions_ >> sampledScalars >> sampledVectors;

        sampledScalars_.setSize(sampledScalars.size());
        forAll(sampledScalars, i)
        {
            sampledScalars_[i].transfer(sampledScalars[i]);
        }
        sampledVectors_.setSize(sampledVectors.size());
        forAll(sampledVectors, i)
        {
            sampledVectors_[i].transfer(sampledVectors[i]);
        }
    }

    // Check state of Istream
    is.check
    (
        "wallBoundedStreamLineParticle::wallBoundedStreamLineParticle"
        "(const Cloud<wallBoundedStreamLineParticle>&, Istream&, bool)"
    );
}


Foam::wallBoundedStreamLineParticle::wallBoundedStreamLineParticle
(
    const wallBoundedStreamLineParticle& p
)
:
    particle(p),
    meshEdgeStart_(p.meshEdgeStart_),
    diagEdge_(p.diagEdge_),
    lifeTime_(p.lifeTime_),
    sampledPositions_(p.sampledPositions_),
    sampledScalars_(p.sampledScalars_),
    sampledVectors_(p.sampledVectors_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::wallBoundedStreamLineParticle::move
(
    trackingData& td,
    const scalar trackTime
)
{
    wallBoundedStreamLineParticle& p = static_cast
    <
        wallBoundedStreamLineParticle&
    >(*this);


    // Check position is inside tet
    //checkInside();

    td.switchProcessor = false;
    td.keepParticle = true;

    scalar tEnd = (1.0 - stepFraction())*trackTime;
    scalar maxDt = mesh_.bounds().mag();

    while
    (
        td.keepParticle
    && !td.switchProcessor
    && lifeTime_ > 0
    )
    {
        // set the lagrangian time-step
        scalar dt = maxDt;

        --lifeTime_;

        // Get sampled velocity and fields. Store if position changed.
        vector U = sample(td);

        // !user parameter!
        if (dt < SMALL)
        {
            // Force removal
            lifeTime_ = 0;
            break;
        }


        if (td.trackLength_ < GREAT)
        {
            dt = td.trackLength_;
        }


        scalar fraction = trackToEdge(td, position() + dt*U);
        dt *= fraction;

        tEnd -= dt;
        stepFraction() = 1.0 - tEnd/trackTime;


        if (tEnd <= ROOTVSMALL)
        {
            // Force removal
            lifeTime_ = 0;
        }

        if
        (
            !td.keepParticle
        ||  td.switchProcessor
        ||  lifeTime_ == 0
        )
        {
            break;
        }
    }


    if (!td.keepParticle || lifeTime_ == 0)
    {
        if (lifeTime_ == 0)
        {
            if (debug)
            {
                Pout<< "wallBoundedStreamLineParticle :"
                    << " Removing stagnant particle:"
                    << p.position()
                    << " sampled positions:" << sampledPositions_.size()
                    << endl;
            }
            td.keepParticle = false;
        }
        else
        {
            // Normal exit. Store last position and fields
            sample(td);

            if (debug)
            {
                Pout<< "wallBoundedStreamLineParticle : Removing particle:"
                    << p.position()
                    << " sampled positions:" << sampledPositions_.size()
                    << endl;
            }
        }

        // Transfer particle data into trackingData.
        {
            //td.allPositions_.append(sampledPositions_);
            td.allPositions_.append(vectorList());
            vectorList& top = td.allPositions_.last();
            top.transfer(sampledPositions_);
        }

        forAll(sampledScalars_, i)
        {
            //td.allScalars_[i].append(sampledScalars_[i]);
            td.allScalars_[i].append(scalarList());
            scalarList& top = td.allScalars_[i].last();
            top.transfer(sampledScalars_[i]);
        }
        forAll(sampledVectors_, i)
        {
            //td.allVectors_[i].append(sampledVectors_[i]);
            td.allVectors_[i].append(vectorList());
            vectorList& top = td.allVectors_[i].last();
            top.transfer(sampledVectors_[i]);
        }
    }

    return td.keepParticle;
}


bool Foam::wallBoundedStreamLineParticle::hitPatch
(
    const polyPatch&,
    trackingData& td,
    const label patchI,
    const scalar trackFraction,
    const tetIndices& tetIs
)
{
    // Disable generic patch interaction
    return false;
}


void Foam::wallBoundedStreamLineParticle::hitWedgePatch
(
    const wedgePolyPatch& pp,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::wallBoundedStreamLineParticle::hitSymmetryPatch
(
    const symmetryPolyPatch& pp,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::wallBoundedStreamLineParticle::hitCyclicPatch
(
    const cyclicPolyPatch& pp,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::wallBoundedStreamLineParticle::hitProcessorPatch
(
    const processorPolyPatch& pp,
    trackingData& td
)
{
    // Switch particle
    td.switchProcessor = true;

    // Adapt edgeStart_ for other side.
    // E.g. if edgeStart_ is 1 then the edge is between vertex 1 and 2 so
    // on the other side between 2 and 3 so edgeStart_ should be
    // f.size()-edgeStart_-1.

    const Foam::face& f = mesh_.faces()[face()];

    if (meshEdgeStart_ != -1)
    {
        meshEdgeStart_ = f.size()-meshEdgeStart_-1;
    }
    else
    {
        // diagEdge_ is relative to faceBasePt
        diagEdge_ = f.size()-diagEdge_;
    }
}


void Foam::wallBoundedStreamLineParticle::hitWallPatch
(
    const wallPolyPatch& wpp,
    trackingData& td,
    const tetIndices&
)
{}


void Foam::wallBoundedStreamLineParticle::hitPatch
(
    const polyPatch& wpp,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::wallBoundedStreamLineParticle::readFields
(
    Cloud<wallBoundedStreamLineParticle>& c
)
{
    if (!c.size())
    {
        return;
    }

    particle::readFields(c);

    IOField<label> meshEdgeStart
    (
        c.fieldIOobject("meshEdgeStart", IOobject::MUST_READ)
    );

    IOField<label> diagEdge
    (
        c.fieldIOobject("diagEdge_", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, diagEdge);

    IOField<label> lifeTime
    (
        c.fieldIOobject("lifeTime", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, lifeTime);

    vectorFieldIOField sampledPositions
    (
        c.fieldIOobject("sampledPositions", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, sampledPositions);

    label i = 0;
    forAllIter(Cloud<wallBoundedStreamLineParticle>, c, iter)
    {
        iter().meshEdgeStart_ = meshEdgeStart[i];
        iter().diagEdge_ = diagEdge[i];
        iter().lifeTime_ = lifeTime[i];
        iter().sampledPositions_.transfer(sampledPositions[i]);
        i++;
    }
}


void Foam::wallBoundedStreamLineParticle::writeFields
(
    const Cloud<wallBoundedStreamLineParticle>& c
)
{
    particle::writeFields(c);

    label np =  c.size();

    IOField<label> meshEdgeStart
    (
        c.fieldIOobject("meshEdgeStart", IOobject::NO_READ),
        np
    );
    IOField<label> diagEdge
    (
        c.fieldIOobject("diagEdge", IOobject::NO_READ),
        np
    );
    IOField<label> lifeTime
    (
        c.fieldIOobject("lifeTime", IOobject::NO_READ),
        np
    );
    vectorFieldIOField sampledPositions
    (
        c.fieldIOobject("sampledPositions", IOobject::NO_READ),
        np
    );

    label i = 0;
    forAllConstIter(Cloud<wallBoundedStreamLineParticle>, c, iter)
    {
        meshEdgeStart[i] = iter().meshEdgeStart_;
        diagEdge[i] = iter().diagEdge_;
        lifeTime[i] = iter().lifeTime_;
        sampledPositions[i] = iter().sampledPositions_;
        i++;
    }

    meshEdgeStart.write();
    diagEdge.write();
    lifeTime.write();
    sampledPositions.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const wallBoundedStreamLineParticle& p
)
{
    os  << static_cast<const particle&>(p)
        << token::SPACE << p.meshEdgeStart_
        << token::SPACE << p.diagEdge_
        << token::SPACE << p.lifeTime_
        << token::SPACE << p.sampledPositions_
        << token::SPACE << p.sampledScalars_
        << token::SPACE << p.sampledVectors_;

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const wallBoundedStreamLineParticle&)"
    );

    return os;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<wallBoundedStreamLineParticle>& ip
)
{
    const wallBoundedStreamLineParticle& p = ip.t_;

    tetPointRef tpr(p.currentTet());

    os  << "    " << static_cast<const particle&>(p) << nl
        << "    tet:" << nl;
    os  << "    ";
    meshTools::writeOBJ(os, tpr.a());
    os  << "    ";
    meshTools::writeOBJ(os, tpr.b());
    os  << "    ";
    meshTools::writeOBJ(os, tpr.c());
    os  << "    ";
    meshTools::writeOBJ(os, tpr.d());
    os  << "    l 1 2" << nl
        << "    l 1 3" << nl
        << "    l 1 4" << nl
        << "    l 2 3" << nl
        << "    l 2 4" << nl
        << "    l 3 4" << nl;
    os  << "    ";
    meshTools::writeOBJ(os, p.position());

    return os;
}



// ************************************************************************* //
