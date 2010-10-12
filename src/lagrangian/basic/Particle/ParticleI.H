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

#include "polyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ParticleType>
inline void Foam::Particle<ParticleType>::findTris
(
    const vector& position,
    DynamicList<label>& faceList,
    const tetPointRef& tet,
    const FixedList<vector, 4>& tetAreas,
    const FixedList<label, 4>& tetPlaneBasePtIs
) const
{
    faceList.clear();

    const point Ct = tet.centre();

    for (label i = 0; i < 4; i++)
    {
        scalar lambda = tetLambda
        (
            Ct,
            position,
            i,
            tetAreas[i],
            tetPlaneBasePtIs[i],
            cellI_,
            tetFaceI_,
            tetPtI_
        );

        if ((lambda > 0.0) && (lambda < 1.0))
        {
            faceList.append(i);
        }
    }
}


template<class ParticleType>
inline Foam::scalar Foam::Particle<ParticleType>::tetLambda
(
    const vector& from,
    const vector& to,
    const label triI,
    const vector& n,
    const label tetPlaneBasePtI,
    const label cellI,
    const label tetFaceI,
    const label tetPtI
) const
{
    const polyMesh& mesh = cloud_.polyMesh_;
    const pointField& pPts = mesh.points();

    if (mesh.moving())
    {
        return movingTetLambda
        (
            from,
            to,
            triI,
            n,
            tetPlaneBasePtI,
            cellI,
            tetFaceI,
            tetPtI
        );
    }

    const point& base = pPts[tetPlaneBasePtI];

    scalar lambdaNumerator = (base - from) & n;
    scalar lambdaDenominator = (to - from) & n;

    if (mag(lambdaDenominator) < SMALL)
    {
        if (mag(lambdaNumerator) < SMALL)
        {
            // Track starts on the face, and is potentially
            // parallel to it.  +-SMALL)/+-SMALL is not a good
            // comparison, return 0.0, in anticipation of tet
            // centre correction.

            return 0.0;
        }
        else
        {
            if (mag((to - from)) < SMALL)
            {
                // Zero length track, not along the face, face
                // cannot be crossed.

                return GREAT;
            }
            else
            {
                // Trajectory is non-zero and parallel to face

                lambdaDenominator = sign(lambdaDenominator)*SMALL;
            }
        }
    }

    return lambdaNumerator/lambdaDenominator;
}



template<class ParticleType>
inline Foam::scalar Foam::Particle<ParticleType>::movingTetLambda
(
    const vector& from,
    const vector& to,
    const label triI,
    const vector& n,
    const label tetPlaneBasePtI,
    const label cellI,
    const label tetFaceI,
    const label tetPtI
) const
{
    const polyMesh& mesh = cloud_.polyMesh_;
    const pointField& pPts = mesh.points();
    const pointField& oldPPts = mesh.oldPoints();

    // Base point of plane at end of motion
    const point& b = pPts[tetPlaneBasePtI];

    // n: Normal of plane at end of motion

    // Base point of plane at start of timestep
    const point& b00 = oldPPts[tetPlaneBasePtI];

    // Base point of plane at start of tracking portion (cast forward by
    // stepFraction)
    point b0 = b00 + stepFraction_*(b - b00);

    // Normal of plane at start of tracking portion
    vector n0 = vector::zero;

    {
        tetIndices tetIs(cellI, tetFaceI, tetPtI, mesh);

        // tet at timestep start
        tetPointRef tet00 = tetIs.oldTet(mesh);

        // tet at timestep end
        tetPointRef tet = tetIs.tet(mesh);

        point tet0PtA = tet00.a() + stepFraction_*(tet.a() - tet00.a());
        point tet0PtB = tet00.b() + stepFraction_*(tet.b() - tet00.b());
        point tet0PtC = tet00.c() + stepFraction_*(tet.c() - tet00.c());
        point tet0PtD = tet00.d() + stepFraction_*(tet.d() - tet00.d());

        // Tracking portion start tet (cast forward by stepFraction)
        tetPointRef tet0(tet0PtA, tet0PtB, tet0PtC, tet0PtD);

        switch (triI)
        {
            case 0:
            {
                n0 = tet0.Sa();
                break;
            }
            case 1:
            {
                n0 = tet0.Sb();
                break;
            }
            case 2:
            {
                n0 = tet0.Sc();
                break;
            }
            case 3:
            {
                n0 = tet0.Sd();
                break;
            }
            default:
            {
                break;
            }
        }
    }

    if (mag(n0) < SMALL)
    {
        // If the old normal is zero (for example in layer addition)
        // then use the current normal;

        n0 = n;
    }

    scalar lambdaNumerator = 0;
    scalar lambdaDenominator = 0;

    vector dP = to - from;
    vector dN = n - n0;
    vector dB = b - b0;
    vector dS = from - b0;

    if (mag(dN) > SMALL)
    {
        scalar a = (dP - dB) & dN;
        scalar b = ((dP - dB) & n0) + (dS & dN);
        scalar c = dS & n0;

        if (mag(a) > SMALL)
        {

            // Solve quadratic for lambda
            scalar discriminant = sqr(b) - 4.0*a*c;

            if (discriminant < 0)
            {
                // Imaginary roots only - face not crossed
                return GREAT;
            }
            else
            {
                // Numerical Recipes in C,
                // Second Edition (1992),
                // Section 5.6.
                // q = -0.5*(b + sgn(b)*sqrt(b^2 - 4ac))
                // x1 = q/a
                // x2 = c/q

                scalar q = -0.5*(b + sign(b)*Foam::sqrt(discriminant));

                if (mag(q) < VSMALL)
                {
                    // If q is zero, then l1 = q/a is the required
                    // value of lambda, and is zero.

                    return 0.0;
                }

                scalar l1 = q/a;
                scalar l2 = c/q;

                // There will be two roots, a big one and a little
                // one, choose the little one.

                if (mag(l1) < mag(l2))
                {
                    return l1;
                }
                else
                {
                    return l2;
                }
            }
        }
        {
            // When a is zero, solve the first order polynomial

            lambdaNumerator = -c;
            lambdaDenominator = b;
        }
    }
    else
    {
        // when n = n0 is zero, there is no plane rotation, solve the
        // first order polynomial

        lambdaNumerator = -(dS & n0);
        lambdaDenominator = ((dP - dB) & n0);

    }

    if (mag(lambdaDenominator) < SMALL)
    {
        if (mag(lambdaNumerator) < SMALL)
        {
            // Track starts on the face, and is potentially
            // parallel to it.  +-SMALL)/+-SMALL is not a good
            // comparison, return 0.0, in anticipation of tet
            // centre correction.

            return 0.0;
        }
        else
        {
            if (mag((to - from)) < SMALL)
            {
                // Zero length track, not along the face, face
                // cannot be crossed.

                return GREAT;
            }
            else
            {
                // Trajectory is non-zero and parallel to face

                lambdaDenominator = sign(lambdaDenominator)*SMALL;
            }
        }
    }

    return lambdaNumerator/lambdaDenominator;
}


template<class ParticleType>
inline void Foam::Particle<ParticleType>::tetNeighbour(label triI)
{
    const polyMesh& mesh = cloud_.polyMesh_;

    const labelList& pOwner = mesh.faceOwner();
    const faceList& pFaces = mesh.faces();

    bool own = (pOwner[tetFaceI_] == cellI_);

    const Foam::face& f = pFaces[tetFaceI_];

    label tetBasePtI = mesh.tetBasePtIs()[tetFaceI_];

    label facePtI = (tetPtI_ + tetBasePtI) % f.size();
    label otherFacePtI = f.fcIndex(facePtI);

    if (tetBasePtI == -1)
    {
        FatalErrorIn
        (
            "inline void Foam::Particle<ParticleType>::tetNeighbour(label triI)"
        )
            << "No base point for face " << tetFaceI_ << ", " << f
            << ", produces a valid tet decomposition."
            << abort(FatalError);
    }

    switch (triI)
    {
        case 0:
        {
            // Crossing this triangle changes tet to that in the
            // neighbour cell over tetFaceI

            // Modification of cellI_ will happen by other indexing,
            // tetFaceI_ and tetPtI don't change.

            break;
        }
        case 1:
        {
            crossEdgeConnectedFace
            (
                cellI_,
                tetFaceI_,
                tetPtI_,
                Foam::edge(f[facePtI], f[otherFacePtI])
            );

            break;
        }
        case 2:
        {
            if (own)
            {
                if (tetPtI_ < f.size() - 2)
                {
                    tetPtI_ = f.fcIndex(tetPtI_);
                }
                else
                {
                    crossEdgeConnectedFace
                    (
                        cellI_,
                        tetFaceI_,
                        tetPtI_,
                        Foam::edge(f[tetBasePtI], f[otherFacePtI])
                    );
                }
            }
            else
            {
                if (tetPtI_ > 1)
                {
                    tetPtI_ = f.rcIndex(tetPtI_);
                }
                else
                {
                    crossEdgeConnectedFace
                    (
                        cellI_,
                        tetFaceI_,
                        tetPtI_,
                        Foam::edge(f[tetBasePtI], f[facePtI])
                    );
                }
            }

            break;
        }
        case 3:
        {
            if (own)
            {
                if (tetPtI_ > 1)
                {
                    tetPtI_ = f.rcIndex(tetPtI_);
                }
                else
                {
                    crossEdgeConnectedFace
                    (
                        cellI_,
                        tetFaceI_,
                        tetPtI_,
                        Foam::edge(f[tetBasePtI], f[facePtI])
                    );
                }
            }
            else
            {
                if (tetPtI_ < f.size() - 2)
                {
                    tetPtI_ = f.fcIndex(tetPtI_);
                }
                else
                {
                    crossEdgeConnectedFace
                    (
                        cellI_,
                        tetFaceI_,
                        tetPtI_,
                        Foam::edge(f[tetBasePtI], f[otherFacePtI])
                    );
                }
            }

            break;
        }
        default:
        {
            FatalErrorIn
            (
                "inline void "
                "Foam::Particle<ParticleType>::tetNeighbour(label triI)"
            )
                << "Tet tri face index error, can only be 0..3, supplied "
                << triI << nl
                << abort(FatalError);

            break;
        }
    }
}

template<class ParticleType>
inline void Foam::Particle<ParticleType>::crossEdgeConnectedFace
(
    const label& cellI,
    label& tetFaceI,
    label& tetPtI,
    const edge& e
)
{
    const polyMesh& mesh = cloud_.polyMesh_;

    const faceList& pFaces = mesh.faces();
    const cellList& pCells = mesh.cells();

    const Foam::face& f = pFaces[tetFaceI];

    const Foam::cell& thisCell = pCells[cellI];

    forAll(thisCell, cFI)
    {
        // Loop over all other faces of this cell and
        // find the one that shares this edge

        label fI = thisCell[cFI];

        if (tetFaceI == fI)
        {
            continue;
        }

        const Foam::face& otherFace = pFaces[fI];

        label edDir = otherFace.edgeDirection(e);

        if (edDir == 0)
        {
            continue;
        }
        else
        {
            //Found edge on other face
            tetFaceI = fI;

            label eIndex = -1;

            if (edDir == 1)
            {
                // Edge is in the forward circulation of this face, so
                // work with the start point of the edge

                eIndex = findIndex(otherFace, e.start());
            }
            else
            {
                // edDir == -1, so the edge is in the reverse
                // circulation of this face, so work with the end
                // point of the edge

                eIndex = findIndex(otherFace, e.end());
            }

            label tetBasePtI = mesh.tetBasePtIs()[fI];

            if (tetBasePtI == -1)
            {
                FatalErrorIn
                (
                    "inline void "
                    "Foam::Particle<ParticleType>::crossEdgeConnectedFace"
                    "("
                        "const label& cellI,"
                        "label& tetFaceI,"
                        "label& tetPtI,"
                        "const edge& e"
                    ")"
                )
                    << "No base point for face " << fI << ", " << f
                    << ", produces a decomposition that has a minimum "
                    << "volume greater than tolerance."
                    << abort(FatalError);
            }

            // Find eIndex relative to the base point on new face
            eIndex -= tetBasePtI;

            if (neg(eIndex))
            {
                eIndex = (eIndex + otherFace.size()) % otherFace.size();
            }

            if (eIndex == 0)
            {
                // The point is the base point, so this is first tet
                // in the face circulation

                tetPtI = 1;
            }
            else if (eIndex == otherFace.size() - 1)
            {
                // The point is the last before the base point, so
                // this is the last tet in the face circulation

                tetPtI = otherFace.size() - 2;
            }
            else
            {
                tetPtI = eIndex;
            }

            break;
        }
    }
}


template<class ParticleType>
inline void Foam::Particle<ParticleType>::hitWallFaces
(
    const vector& from,
    const vector& to,
    scalar& lambdaMin,
    tetIndices& closestTetIs
)
{
    if (!(cloud_.hasWallImpactDistance() && cloud_.cellHasWallFaces()[cellI_]))
    {
        return;
    }

    const polyMesh& mesh = cloud_.polyMesh_;
    const faceList& pFaces = mesh.faces();

    const ParticleType& p = static_cast<const ParticleType&>(*this);

    const Foam::cell& thisCell = mesh.cells()[cellI_];

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(thisCell, cFI)
    {
        label fI = thisCell[cFI];

        if (cloud_.internalFace(fI))
        {
            continue;
        }

        label patchI = patches.patchID()[fI - mesh.nInternalFaces()];

        if (isA<wallPolyPatch>(patches[patchI]))
        {
            // Get the decomposition of this wall face

            const List<tetIndices>& faceTetIs =
                polyMeshTetDecomposition::faceTetIndices(mesh, fI, cellI_);

            const Foam::face& f = pFaces[fI];

            forAll(faceTetIs, tI)
            {
                const tetIndices& tetIs = faceTetIs[tI];

                triPointRef tri = tetIs.faceTri(mesh);

                vector n = tri.normal();

                vector nHat = n/mag(n);

                // Radius of particle with respect to this wall face
                // triangle.  Assuming that the wallImpactDistance
                // does not change as the particle or the mesh moves
                // forward in time.
                scalar r = p.wallImpactDistance(nHat);

                vector toPlusRNHat = to + r*nHat;

                // triI = 0 because it is the cell face tri of the tet
                // we are concerned with.
                scalar tetClambda = tetLambda
                (
                    tetIs.tet(mesh).centre(),
                    toPlusRNHat,
                    0,
                    n,
                    f[tetIs.faceBasePt()],
                    cellI_,
                    fI,
                    tetIs.tetPt()
                );

                if ((tetClambda <= 0.0) || (tetClambda >= 1.0))
                {
                    // toPlusRNHat is not on the outside of the plane of
                    // the wall face tri, the tri cannot be hit.
                    continue;
                }

                // Check if the actual trajectory of the near-tri
                // points intersects the triangle.

                vector fromPlusRNHat = from + r*nHat;

                // triI = 0 because it is the cell face tri of the tet
                // we are concerned with.
                scalar lambda = tetLambda
                (
                    fromPlusRNHat,
                    toPlusRNHat,
                    0,
                    n,
                    f[tetIs.faceBasePt()],
                    cellI_,
                    fI,
                    tetIs.tetPt()
                );

                pointHit hitInfo(vector::zero);

                if (mesh.moving())
                {
                    // For a moving mesh, the position of wall
                    // triangle needs to be moved in time to be
                    // consistent with the moment defined by the
                    // current value of stepFraction and the value of
                    // lambda just calculated.

                    // Total fraction thought the timestep of the
                    // motion, including stepFraction before the
                    // current tracking step and the current
                    // lambda
                    // i.e.
                    // let s = stepFraction, l = lambda
                    // Motion of x in time:
                    // |-----------------|---------|---------|
                    // x00               x0        xi        x
                    //
                    // where xi is the correct value of x at the required
                    // tracking instant.
                    //
                    // x0 = x00 + s*(x - x00) = s*x + (1 - s)*x00
                    //
                    // i.e. the motion covered by previous tracking portions
                    // within this timestep, and
                    //
                    // xi = x0 + l*(x - x0)
                    //    = l*x + (1 - l)*x0
                    //    = l*x + (1 - l)*(s*x + (1 - s)*x00)
                    //    = (s + l - s*l)*x + (1 - (s + l - s*l))*x00
                    //
                    // let m = (s + l - s*l)
                    //
                    // xi = m*x + (1 - m)*x00 = x00 + m*(x - x00);
                    //
                    // In the same form as before.

                    // Clip lambda to 0.0-1.0 to ensure that sensible
                    // positions are used for triangle intersections.
                    scalar lam = max(0.0, min(1.0, lambda));

                    scalar m = stepFraction_ + lam - (stepFraction_*lam);

                    triPointRef tri00 = tetIs.oldFaceTri(mesh);

                    // Use SMALL positive tolerance to make the triangle
                    // slightly "fat" to improve robustness.  Intersection
                    // is calculated as the ray (from + r*nHat) -> (to +
                    // r*nHat).

                    point tPtA = tri00.a() + m*(tri.a() - tri00.a());
                    point tPtB = tri00.b() + m*(tri.b() - tri00.b());
                    point tPtC = tri00.c() + m*(tri.c() - tri00.c());

                    triPointRef t(tPtA, tPtB, tPtC);

                    // The point fromPlusRNHat + m*(to - from) is on the
                    // plane of the triangle.  Determine the
                    // intersection with this triangle by testing if
                    // this point is inside or outside of the triangle.
                    hitInfo = t.intersection
                    (
                        fromPlusRNHat + m*(to - from),
                        t.normal(),
                        intersection::FULL_RAY,
                        SMALL
                    );
                }
                else
                {
                    // Use SMALL positive tolerance to make the triangle
                    // slightly "fat" to improve robustness.  Intersection
                    // is calculated as the ray (from + r*nHat) -> (to +
                    // r*nHat).
                    hitInfo = tri.intersection
                    (
                        fromPlusRNHat,
                        (to - from),
                        intersection::FULL_RAY,
                        SMALL
                    );
                }

                if (hitInfo.hit())
                {
                    if (lambda < lambdaMin)
                    {
                        lambdaMin = lambda;

                        faceI_ = fI;

                        closestTetIs = tetIs;
                    }
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
inline Foam::Particle<ParticleType>::trackData::trackData
(
    Cloud<ParticleType>& cloud
)
:
    cloud_(cloud)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
inline Foam::Cloud<ParticleType>&
Foam::Particle<ParticleType>::trackData::cloud()
{
    return cloud_;
}


template<class ParticleType>
inline const Foam::Cloud<ParticleType>&
Foam::Particle<ParticleType>::cloud() const
{
    return cloud_;
}


template<class ParticleType>
inline const Foam::vector& Foam::Particle<ParticleType>::position() const
{
    return position_;
}


template<class ParticleType>
inline Foam::vector& Foam::Particle<ParticleType>::position()
{
    return position_;
}


template<class ParticleType>
inline Foam::label Foam::Particle<ParticleType>::cell() const
{
    return cellI_;
}


template<class ParticleType>
inline Foam::label& Foam::Particle<ParticleType>::cell()
{
    return cellI_;
}


template<class ParticleType>
inline Foam::label Foam::Particle<ParticleType>::tetFace() const
{
    return tetFaceI_;
}


template<class ParticleType>
inline Foam::label& Foam::Particle<ParticleType>::tetFace()
{
    return tetFaceI_;
}


template<class ParticleType>
inline Foam::label Foam::Particle<ParticleType>::tetPt() const
{
    return tetPtI_;
}


template<class ParticleType>
inline Foam::label& Foam::Particle<ParticleType>::tetPt()
{
    return tetPtI_;
}


template<class ParticleType>
inline Foam::tetIndices Foam::Particle<ParticleType>::currentTetIndices() const
{
    return tetIndices(cellI_, tetFaceI_, tetPtI_, cloud_.polyMesh_);
}


template<class ParticleType>
inline Foam::tetPointRef Foam::Particle<ParticleType>::currentTet() const
{
    return currentTetIndices().tet(cloud_.polyMesh_);
}


template<class ParticleType>
inline Foam::vector Foam::Particle<ParticleType>::normal() const
{
    return currentTetIndices().faceTri(cloud_.polyMesh_).normal();
}


template<class ParticleType>
inline Foam::vector
Foam::Particle<ParticleType>::oldNormal() const
{
    return currentTetIndices().oldFaceTri(cloud_.polyMesh_).normal();
}


template<class ParticleType>
inline Foam::label Foam::Particle<ParticleType>::face() const
{
    return faceI_;
}


template<class ParticleType>
inline Foam::label& Foam::Particle<ParticleType>::face()
{
    return faceI_;
}


template<class ParticleType>
inline void Foam::Particle<ParticleType>::initCellFacePt()
{
    if (cellI_ == -1)
    {
        cloud_.findCellFacePt
        (
            position_,
            cellI_,
            tetFaceI_,
            tetPtI_
        );

        if (cellI_ == -1)
        {
            FatalErrorIn
            (
                "Foam::Particle<ParticleType>::Particle"
                "("
                    "void Foam::Particle<ParticleType>::initCellFacePt()"
                ")"
            )   << "cell, tetFace and tetPt search failure at position "
                << position_ << nl
                << abort(FatalError);
        }
    }
    else
    {
        cloud_.findFacePt(cellI_, position_, tetFaceI_, tetPtI_);

        if (tetFaceI_ == -1 || tetPtI_ == -1)
        {
            label oldCellI = cellI_;

            cloud_.findCellFacePt
            (
                position_,
                cellI_,
                tetFaceI_,
                tetPtI_
            );

            if (cellI_ == -1 || tetFaceI_ == -1 || tetPtI_ == -1)
            {
                // The particle has entered this function with a cell
                // number, but hasn't been able to find a cell to
                // occupy.

                if(!cloud_.polyMesh_.pointInCellBB(position_, oldCellI, 0.1))
                {
                    // If the position is not inside the (slightly
                    // extended) bound-box of the cell that it thought
                    // it should be in, then this is considered an
                    // error.

                    FatalErrorIn
                    (
                        "Foam::Particle<ParticleType>::Particle"
                        "("
                            "void "
                            "Foam::Particle<ParticleType>::initCellFacePt()"
                        ")"
                    )   << "cell, tetFace and tetPt search failure at position "
                        << position_ << nl
                        << "for requested cell " << oldCellI << nl
                        << abort(FatalError);
                }

                // The position is in the (slightly extended)
                // bound-box of the cell.  This situation may arise
                // because the face decomposition of the cell is not
                // the same as when the particle acquired the cell
                // index.  For example, it has been read into a mesh
                // that has made a different face base-point decision
                // for a boundary face and now this particle is in a
                // position that is not in the mesh.  Gradually move
                // the particle towards the centre of the cell that it
                // thought that it was in.

                cellI_ = oldCellI;

                point newPosition = position_;

                const point& cC = cloud_.polyMesh_.cellCentres()[cellI_];

                label trap(1.0/Cloud<ParticleType>::trackingCorrectionTol + 1);

                label iterNo = 0;

                do
                {
                    newPosition +=
                        Cloud<ParticleType>::trackingCorrectionTol
                       *(cC - position_);

                    cloud_.findFacePt
                    (
                        cellI_,
                        newPosition,
                        tetFaceI_,
                        tetPtI_
                    );

                    iterNo++;

                } while (tetFaceI_ < 0  && iterNo <= trap);

                if(tetFaceI_ == -1)
                {
                    FatalErrorIn
                    (
                        "Foam::Particle<ParticleType>::Particle"
                        "("
                        "void "
                        "Foam::Particle<ParticleType>::initCellFacePt()"
                        ")"
                    )   << "cell, tetFace and tetPt search failure at position "
                        << position_ << nl
                        << abort(FatalError);
                }

                WarningIn
                (
                    "Foam::Particle<ParticleType>::Particle"
                    "("
                        "void Foam::Particle<ParticleType>::initCellFacePt()"
                    ")"
                )
                    << "Particle moved from " << position_
                    << " to " << newPosition
                    << " in cell " << cellI_
                    << " tetFace " << tetFaceI_
                    << " tetPt " << tetPtI_ << nl
                    << "    (A fraction of "
                    << 1.0 - mag(cC - newPosition)/mag(cC - position_)
                    << " of the distance to the cell centre)"
                    << " because a decomposition tetFace and tetPt "
                    << "could not be found."
                    << endl;

                position_ = newPosition;
            }

            if (cellI_ != oldCellI)
            {
                WarningIn
                (
                    "Foam::Particle<ParticleType>::Particle"
                    "("
                        "void Foam::Particle<ParticleType>::initCellFacePt()"
                    ")"
                )
                    << "Particle at position " << position_
                    << " searched for a cell, tetFace and tetPt." << nl
                    << "    Found"
                    << " cell " << cellI_
                    << " tetFace " << tetFaceI_
                    << " tetPt " << tetPtI_ << nl
                    << "    This is a different cell to that which was supplied"
                    << " (" << oldCellI << ")." << nl
                    << endl;
            }
        }
    }
}


template<class ParticleType>
inline bool Foam::Particle<ParticleType>::onBoundary() const
{
    return faceI_ != -1 && faceI_ >= cloud_.pMesh().nInternalFaces();
}


template<class ParticleType>
inline Foam::scalar& Foam::Particle<ParticleType>::stepFraction()
{
    return stepFraction_;
}


template<class ParticleType>
inline Foam::scalar Foam::Particle<ParticleType>::stepFraction() const
{
    return stepFraction_;
}


template<class ParticleType>
inline Foam::label Foam::Particle<ParticleType>::origProc() const
{
    return origProc_;
}


template<class ParticleType>
inline Foam::label Foam::Particle<ParticleType>::origId() const
{
    return origId_;
}


template<class ParticleType>
inline bool Foam::Particle<ParticleType>::softImpact() const
{
    return false;
}


template<class ParticleType>
inline Foam::scalar Foam::Particle<ParticleType>::currentTime() const
{
    return
        cloud_.pMesh().time().value()
      + stepFraction_*cloud_.pMesh().time().deltaTValue();
}


template<class ParticleType>
inline Foam::label Foam::Particle<ParticleType>::patch(const label faceI) const
{
    return cloud_.facePatch(faceI);
}


template<class ParticleType>
inline Foam::label Foam::Particle<ParticleType>::patchFace
(
    const label patchI,
    const label faceI
) const
{
    return cloud_.patchFace(patchI, faceI);
}


template<class ParticleType>
inline Foam::scalar
Foam::Particle<ParticleType>::wallImpactDistance(const vector&) const
{
    return 0.0;
}


template<class ParticleType>
inline Foam::label Foam::Particle<ParticleType>::faceInterpolation() const
{
    return faceI_;
}


// ************************************************************************* //
