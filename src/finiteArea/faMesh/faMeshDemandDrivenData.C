/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "faMesh.H"
#include "faMeshLduAddressing.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "fac.H"
#include "cartesianCS.H"
#include "scalarMatrices.H"
#include "processorFaPatch.H"
#include "processorFaPatchFields.H"
#include "emptyFaPatchFields.H"
#include "wedgeFaPatch.H"
#include "triPointRef.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Define an area-weighted normal for three points (defining a triangle)
// (p0, p1, p2) are the base, first, second points respectively
//
// From the original Tukovic code:
//
//   vector n = (d1^d2)/mag(d1^d2);
//   scalar sinAlpha = mag(d1^d2)/(mag(d1)*mag(d2));
//   scalar w = sinAlpha/(mag(d1)*mag(d2));
//
// ie, normal weighted by area, sine angle and inverse distance squared.
// - area : larger weight for larger areas
// - sin  : lower weight for narrow angles (eg, shards)
// - inv distance squared : lower weights for distant points
//
// The above refactored, with 0.5 for area:
//
//   (d1 ^ d2) / (2 * magSqr(d1) * magSqr(d2))

static inline vector areaInvDistSqrWeightedNormal
(
    const vector& a,
    const vector& b
)
{
    const scalar s(2*magSqr(a)*magSqr(b));

    return s < VSMALL ? Zero : (a ^ b) / s;
}


// The area normal for the face dual (around base-point)
// described by the right/left edge points and the centre point
//
// The adjustment for 1/2 edge point (the dual point) is done internally
static inline vector areaInvDistSqrWeightedNormalDualEdge
(
    const point& basePoint,
    const point& rightPoint,
    const point& leftPoint,
    const point& centrePoint
)
{
    const vector mid(centrePoint - basePoint);

    return
    (
        areaInvDistSqrWeightedNormal
        (
            0.5*(rightPoint - basePoint),   // vector to right 1/2 edge
            mid
        )
      + areaInvDistSqrWeightedNormal
        (
            mid,
            0.5*(leftPoint - basePoint)     // vector to left 1/2 edge
        )
    );
}


// Calculate transform tensor with reference vector (unitAxis1)
// and direction vector (axis2).
//
// This is nearly identical to the meshTools axesRotation
// with an E3_E1 transformation with the following exceptions:
//
// - axis1 (e3 == unitAxis1): is already normalized (unit vector)
// - axis2 (e1 == dirn): no difference
// - transformation is row-vectors, not column-vectors
static inline tensor rotation_e3e1
(
    const vector& unitAxis1,
    vector dirn
)
{
    dirn.removeCollinear(unitAxis1);
    dirn.normalise();

    // Set row vectors
    return tensor
    (
        dirn,
        (unitAxis1^dirn),
        unitAxis1
    );
}


// Simple area-weighted normal calculation for the specified edge vector
// and its owner/neighbour face centres (internal edges).
//
// Uses four triangles since adjacent faces may be non-planar
// and/or the edge centre is skewed from the face centres.

/*---------------------------------------------------------------------------*\
 *          *           |
 *         /|\          | Triangles (0,1) are on the owner side.
 *        / | \         |
 *       /  |  \        | Triangles (2,3) are on the neighbour side.
 *      /  1|3  \       |
 * own *----|----* nbr  | Boundary edges are the same, but without a neighbour
 *      \  0|2  /       |
 *       \  |  /        |
 *        \ | /         |
 *         \|/          |
 *          *           |
\*---------------------------------------------------------------------------*/

static inline vector calcEdgeNormalFromFace
(
    const linePointRef& edgeVec,
    const point& ownCentre,
    const point& neiCentre
)
{
    const scalar magEdge(edgeVec.mag());

    if (magEdge < ROOTVSMALL)
    {
        return Zero;
    }

    const vector edgeCtr = edgeVec.centre();

    vector edgeNorm
    (
        // owner
        triPointRef(edgeCtr, ownCentre, edgeVec.first()).areaNormal()
      + triPointRef(edgeCtr, edgeVec.last(), ownCentre).areaNormal()
        // neighbour
      + triPointRef(edgeCtr, edgeVec.first(), neiCentre).areaNormal()
      + triPointRef(edgeCtr, neiCentre, edgeVec.last()).areaNormal()
    );

    // Requires a unit-vector (already tested its mag)
    edgeNorm.removeCollinear(edgeVec.vec()/magEdge);
    edgeNorm.normalise();

    // Scale with the original edge length
    return magEdge*edgeNorm;
}


// As above, but for boundary edgess (no neighbour)
static inline vector calcEdgeNormalFromFace
(
    const linePointRef& edgeVec,
    const point& ownCentre
)
{
    const scalar magEdge(edgeVec.mag());

    if (magEdge < ROOTVSMALL)
    {
        return Zero;
    }

    const vector edgeCtr = edgeVec.centre();

    vector edgeNorm
    (
        // owner
        triPointRef(edgeCtr, ownCentre, edgeVec.first()).areaNormal()
      + triPointRef(edgeCtr, edgeVec.last(), ownCentre).areaNormal()
    );

    // Requires a unit-vector (already tested its mag)
    edgeNorm.removeCollinear(edgeVec.vec()/magEdge);
    edgeNorm.normalise();

    // Scale with the original edge length
    return magEdge*edgeNorm;
}

} // End namespace Foam


namespace Foam
{

// A bitSet (size patch nPoints()) with boundary points marked
static Foam::bitSet markupBoundaryPoints(const uindirectPrimitivePatch& p)
{
    // Initially all unmarked
    bitSet markPoints(p.nPoints());
    for (const edge& e : p.boundaryEdges())
    {
        // Mark boundary points
        markPoints.set(e.first());
        markPoints.set(e.second());
    }

    return markPoints;
}


} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faMesh::calcLduAddressing() const
{
    DebugInFunction
        << "Calculating addressing" << endl;

    if (lduPtr_)
    {
        FatalErrorInFunction
            << "lduPtr_ already allocated"
            << abort(FatalError);
    }

    lduPtr_ = new faMeshLduAddressing(*this);
}


void Foam::faMesh::calcPatchStarts() const
{
    DebugInFunction
        << "Calculating patch starts" << endl;

    if (patchStartsPtr_)
    {
        FatalErrorInFunction
            << "patchStartsPtr_ already allocated"
            << abort(FatalError);
    }

    patchStartsPtr_ = new labelList(boundary().patchStarts());
}


void Foam::faMesh::calcLe() const
{
    DebugInFunction
        << "Calculating local edges" << endl;

    if (LePtr_)
    {
        FatalErrorInFunction
            << "LePtr_ already allocated"
            << abort(FatalError);
    }

    LePtr_ =
        new edgeVectorField
        (
            IOobject
            (
                "Le",
                mesh().pointsInstance(),
                meshSubDir,
                mesh()
            ),
            *this,
            dimLength
        );

    edgeVectorField& Le = *LePtr_;

    const pointField& localPoints = points();

    if (faMesh::geometryOrder() < 1)
    {
        // Simple (primitive) edge normal calculations.
        // These are primarly designed to avoid any communication
        // but are thus necessarily inconsistent across processor boundaries!

        // Reasonable to assume that the volume mesh already has faceCentres
        // eg, used magSf somewhere.
        // Can use these instead of triggering our calcAreaCentres().

        WarningInFunction
            << "Using geometryOrder < 1 : "
               "simplified edge area-normals for Le() calculation"
            << endl;

        UIndirectList<vector> fCentres(mesh().faceCentres(), faceLabels());

        // Flat addressing
        vectorField edgeNormals(nEdges_);

        // Internal (edge normals)
        for (label edgei = 0; edgei < nInternalEdges_; ++edgei)
        {
            edgeNormals[edgei] =
                calcEdgeNormalFromFace
                (
                    edges_[edgei].line(localPoints),
                    fCentres[owner()[edgei]],
                    fCentres[neighbour()[edgei]]
                );
        }

        // Boundary (edge normals). Like above, but only has owner
        for (label edgei = nInternalEdges_; edgei < nEdges_; ++edgei)
        {
            edgeNormals[edgei] =
                calcEdgeNormalFromFace
                (
                    edges_[edgei].line(localPoints),
                    fCentres[owner()[edgei]]
                );
        }


        // Now use these edge normals for calculating Le

        // Internal (edge vector)
        {
            vectorField& internalFld = Le.ref();
            for (label edgei = 0; edgei < nInternalEdges_; ++edgei)
            {
                vector& leVector = internalFld[edgei];

                vector edgeVec = edges_[edgei].vec(localPoints);
                const scalar magEdge(mag(edgeVec));

                if (magEdge < ROOTVSMALL)
                {
                    // Too small
                    leVector = Zero;
                    continue;
                }

                const vector edgeCtr = edges_[edgei].centre(localPoints);
                const vector& edgeNorm = edgeNormals[edgei];
                const vector& ownCentre = fCentres[owner()[edgei]];

                leVector = magEdge*normalised(edgeVec ^ edgeNorm);
                leVector *= sign(leVector & (edgeCtr - ownCentre));
            }
        }

        // Boundary (edge vector)

        forAll(boundary(), patchi)
        {
            const labelUList& bndEdgeFaces = boundary()[patchi].edgeFaces();

            const edgeList::subList bndEdges =
                boundary()[patchi].patchSlice(edges_);

            vectorField& patchLe = Le.boundaryFieldRef()[patchi];

            forAll(patchLe, bndEdgei)
            {
                vector& leVector = patchLe[bndEdgei];
                const label meshEdgei(boundary()[patchi].start() + bndEdgei);

                vector edgeVec = bndEdges[bndEdgei].vec(localPoints);
                const scalar magEdge(mag(edgeVec));

                if (magEdge < ROOTVSMALL)
                {
                    // Too small
                    leVector = Zero;
                    continue;
                }

                const vector edgeCtr = bndEdges[bndEdgei].centre(localPoints);
                const vector& edgeNorm = edgeNormals[meshEdgei];
                const vector& ownCentre = fCentres[bndEdgeFaces[bndEdgei]];

                leVector = magEdge*normalised(edgeVec ^ edgeNorm);
                leVector *= sign(leVector & (edgeCtr - ownCentre));
            }
        }

        // Done
        return;
    }


    // Longer forms.
    // Using edgeAreaNormals, which uses pointAreaNormals (communication!)

    const edgeVectorField& eCentres = edgeCentres();
    const areaVectorField& fCentres = areaCentres();
    const edgeVectorField& edgeNormals = edgeAreaNormals();

    // Internal (edge vector)

    {
        vectorField& internalFld = Le.ref();
        for (label edgei = 0; edgei < nInternalEdges_; ++edgei)
        {
            vector& leVector = internalFld[edgei];

            vector edgeVec = edges_[edgei].vec(localPoints);
            const scalar magEdge(mag(edgeVec));

            if (magEdge < ROOTVSMALL)
            {
                // Too small
                leVector = Zero;
                continue;
            }

            const vector& edgeCtr = eCentres[edgei];
            const vector& edgeNorm = edgeNormals[edgei];
            const vector& ownCentre = fCentres[owner()[edgei]];

            leVector = magEdge*normalised(edgeVec ^ edgeNorm);
            leVector *= sign(leVector & (edgeCtr - ownCentre));
        }
    }

    forAll(boundary(), patchi)
    {
        const labelUList& bndEdgeFaces = boundary()[patchi].edgeFaces();

        const edgeList::subList bndEdges =
            boundary()[patchi].patchSlice(edges_);

        const vectorField& bndEdgeNormals =
            edgeNormals.boundaryField()[patchi];

        vectorField& patchLe = Le.boundaryFieldRef()[patchi];
        const vectorField& patchECentres = eCentres.boundaryField()[patchi];

        forAll(patchLe, bndEdgei)
        {
            vector& leVector = patchLe[bndEdgei];

            vector edgeVec = bndEdges[bndEdgei].vec(localPoints);
            const scalar magEdge(mag(edgeVec));

            if (magEdge < ROOTVSMALL)
            {
                // Too small
                leVector = Zero;
                continue;
            }

            const vector& edgeCtr = patchECentres[bndEdgei];
            const vector& edgeNorm = bndEdgeNormals[bndEdgei];
            const vector& ownCentre = fCentres[bndEdgeFaces[bndEdgei]];

            leVector = magEdge*normalised(edgeVec ^ edgeNorm);
            leVector *= sign(leVector & (edgeCtr - ownCentre));
        }
    }
}


void Foam::faMesh::calcMagLe() const
{
    DebugInFunction
        << "Calculating local edge magnitudes" << endl;

    if (magLePtr_)
    {
        FatalErrorInFunction
            << "magLePtr_ already allocated"
            << abort(FatalError);
    }

    magLePtr_ =
        new edgeScalarField
        (
            IOobject
            (
                "magLe",
                mesh().pointsInstance(),
                meshSubDir,
                mesh()
            ),
            *this,
            dimLength
        );

    edgeScalarField& magLe = *magLePtr_;

    const pointField& localPoints = points();

    // Internal (edge length)
    {
        auto iter = magLe.primitiveFieldRef().begin();

        for (const edge& e : internalEdges())
        {
            *iter = e.mag(localPoints);
            ++iter;
        }
    }

    // Boundary (edge length)
    {
        auto& bfld = magLe.boundaryFieldRef();

        forAll(boundary(), patchi)
        {
            auto iter = bfld[patchi].begin();

            for (const edge& e : boundary()[patchi].patchSlice(edges_))
            {
                *iter = e.mag(localPoints);
                ++iter;
            }
        }
    }
}


void Foam::faMesh::calcAreaCentres() const
{
    DebugInFunction
        << "Calculating face centres" << endl;

    if (centresPtr_)
    {
        FatalErrorInFunction
            << "centresPtr_ already allocated"
            << abort(FatalError);
    }

    centresPtr_ =
        new areaVectorField
        (
            IOobject
            (
                "centres",
                mesh().pointsInstance(),
                meshSubDir,
                mesh()
            ),
            *this,
            dimLength
        );

    areaVectorField& centres = *centresPtr_;

    const pointField& localPoints = points();
    const faceList& localFaces = faces();

    // Internal (face centres)
    // Could also obtain from volume mesh faceCentres()
    forAll(localFaces, facei)
    {
        centres.ref()[facei] = localFaces[facei].centre(localPoints);
    }

    // Boundary (edge centres)
    {
        auto& bfld = centres.boundaryFieldRef();

        forAll(boundary(), patchi)
        {
            auto iter = bfld[patchi].begin();

            for (const edge& e : boundary()[patchi].patchSlice(edges_))
            {
                *iter = e.centre(localPoints);
                ++iter;
            }
        }
    }
}


void Foam::faMesh::calcEdgeCentres() const
{
    DebugInFunction
        << "Calculating edge centres" << endl;

    if (edgeCentresPtr_)
    {
        FatalErrorInFunction
            << "edgeCentresPtr_ already allocated"
            << abort(FatalError);
    }

    edgeCentresPtr_ =
        new edgeVectorField
        (
            IOobject
            (
                "edgeCentres",
                mesh().pointsInstance(),
                meshSubDir,
                mesh()
            ),
            *this,
            dimLength
        );

    edgeVectorField& centres = *edgeCentresPtr_;

    const pointField& localPoints = points();

    // Internal (edge centres)
    {
        auto iter = centres.primitiveFieldRef().begin();

        for (const edge& e : internalEdges())
        {
            *iter = e.centre(localPoints);
            ++iter;
        }
    }

    // Boundary (edge centres)
    {
        auto& bfld = centres.boundaryFieldRef();

        forAll(boundary(), patchi)
        {
            auto iter = bfld[patchi].begin();

            for (const edge& e : boundary()[patchi].patchSlice(edges_))
            {
                *iter = e.centre(localPoints);
                ++iter;
            }
        }
    }
}


void Foam::faMesh::calcS() const
{
    DebugInFunction
        << "Calculating areas" << endl;

    if (SPtr_)
    {
        FatalErrorInFunction
            << "SPtr_ already allocated"
            << abort(FatalError);
    }

    SPtr_ = new DimensionedField<scalar, areaMesh>
    (
        IOobject
        (
            "S",
            time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        dimArea
    );
    auto& S = *SPtr_;

    const pointField& localPoints = points();
    const faceList& localFaces = faces();

    // Could also obtain from volume mesh faceAreas()
    forAll(S, facei)
    {
        S[facei] = localFaces[facei].mag(localPoints);
    }
}


void Foam::faMesh::calcFaceAreaNormals() const
{
    DebugInFunction
         << "Calculating face area normals" << endl;

    if (faceAreaNormalsPtr_)
    {
        FatalErrorInFunction
            << "faceAreaNormalsPtr_ already allocated"
            << abort(FatalError);
    }

    faceAreaNormalsPtr_ =
        new areaVectorField
        (
            IOobject
            (
                "faceAreaNormals",
                mesh().pointsInstance(),
                meshSubDir,
                mesh()
            ),
            *this,
            dimless
        );

    areaVectorField& faceNormals = *faceAreaNormalsPtr_;

    const pointField& localPoints = points();
    const faceList& localFaces = faces();

    // Internal (faces)
    // Could also obtain from volume mesh Sf() + normalise
    vectorField& nInternal = faceNormals.ref();
    forAll(localFaces, faceI)
    {
        nInternal[faceI] = localFaces[faceI].unitNormal(localPoints);
    }

    // Boundary - copy from edges

    const auto& edgeNormalsBoundary = edgeAreaNormals().boundaryField();

    forAll(boundary(), patchI)
    {
        faceNormals.boundaryFieldRef()[patchI] = edgeNormalsBoundary[patchI];
    }
}


void Foam::faMesh::calcEdgeAreaNormals() const
{
    DebugInFunction
        << "Calculating edge area normals" << endl;

    if (edgeAreaNormalsPtr_)
    {
        FatalErrorInFunction
            << "edgeAreaNormalsPtr_ already allocated"
            << abort(FatalError);
    }

    edgeAreaNormalsPtr_ =
        new edgeVectorField
        (
            IOobject
            (
                "edgeAreaNormals",
                mesh().pointsInstance(),
                meshSubDir,
                mesh()
            ),
            *this,
            dimless
        );
    edgeVectorField& edgeAreaNormals = *edgeAreaNormalsPtr_;


    // Starting from point area normals
    const vectorField& pointNormals = pointAreaNormals();


    // Internal edges
    forAll(edgeAreaNormals.internalField(), edgei)
    {
        const edge& e = edges_[edgei];
        const vector edgeVec = e.unitVec(points());

        vector& edgeNorm = edgeAreaNormals.ref()[edgei];

        // Average of both ends
        edgeNorm = (pointNormals[e.first()] + pointNormals[e.second()]);

        edgeNorm.removeCollinear(edgeVec);
        edgeNorm.normalise();
    }

    // Boundary
    auto& bfld = edgeAreaNormals.boundaryFieldRef();

    forAll(boundary(), patchi)
    {
        auto& pfld = bfld[patchi];

        const edgeList::subList patchEdges =
            boundary()[patchi].patchSlice(edges_);

        forAll(patchEdges, bndEdgei)
        {
            const edge& e = patchEdges[bndEdgei];
            const vector edgeVec = e.unitVec(points());

            vector& edgeNorm = pfld[bndEdgei];

            // Average of both ends
            edgeNorm = (pointNormals[e.first()] + pointNormals[e.second()]);

            edgeNorm.removeCollinear(edgeVec);
            edgeNorm.normalise();
        }
    }
}


void Foam::faMesh::calcFaceCurvatures() const
{
    DebugInFunction
        << "Calculating face curvatures" << endl;

    if (faceCurvaturesPtr_)
    {
        FatalErrorInFunction
            << "faceCurvaturesPtr_ already allocated"
            << abort(FatalError);
    }

    faceCurvaturesPtr_ =
        new areaScalarField
        (
            IOobject
            (
                "faceCurvatures",
                mesh().pointsInstance(),
                meshSubDir,
                mesh()
            ),
            *this,
            dimless/dimLength
        );

    areaScalarField& faceCurvatures = *faceCurvaturesPtr_;


//     faceCurvatures =
//         fac::edgeIntegrate(Le()*edgeLengthCorrection())
//         &faceAreaNormals();

    areaVectorField kN(fac::edgeIntegrate(Le()*edgeLengthCorrection()));

    faceCurvatures = sign(kN&faceAreaNormals())*mag(kN);
}


void Foam::faMesh::calcEdgeTransformTensors() const
{
    DebugInFunction
        << "Calculating edge transformation tensors" << endl;

    if (edgeTransformTensorsPtr_)
    {
        FatalErrorInFunction
            << "edgeTransformTensorsPtr_ already allocated"
            << abort(FatalError);
    }

    edgeTransformTensorsPtr_ = new FieldField<Field, tensor>(nEdges_);
    auto& edgeTransformTensors = *edgeTransformTensorsPtr_;

    // Initialize all transformation tensors
    for (label edgei = 0; edgei < nEdges_; ++edgei)
    {
        edgeTransformTensors.set(edgei, new Field<tensor>(3, tensor::I));
    }

    const areaVectorField& Nf = faceAreaNormals();
    const areaVectorField& Cf = areaCentres();

    const edgeVectorField& Ne = edgeAreaNormals();
    const edgeVectorField& Ce = edgeCentres();

    // Internal edges transformation tensors
    for (label edgei = 0; edgei < nInternalEdges_; ++edgei)
    {
        const label ownFacei = owner()[edgei];
        const label neiFacei = neighbour()[edgei];
        auto& tensors = edgeTransformTensors[edgei];

        vector edgeCtr = Ce.internalField()[edgei];

        if (skew())
        {
            edgeCtr -= skewCorrectionVectors().internalField()[edgei];
        }

        // Edge transformation tensor
        tensors[0] =
            rotation_e3e1
            (
                Ne.internalField()[edgei],
                (edgeCtr - Cf.internalField()[ownFacei])
            );

        // Owner transformation tensor
        tensors[1] =
            rotation_e3e1
            (
                Nf.internalField()[ownFacei],
                (edgeCtr - Cf.internalField()[ownFacei])
            );

        // Neighbour transformation tensor
        tensors[2] =
            rotation_e3e1
            (
                Nf.internalField()[neiFacei],
                (Cf.internalField()[neiFacei] - edgeCtr)
            );
    }

    // Boundary edges transformation tensors
    forAll(boundary(), patchi)
    {
        const labelUList& edgeFaces = boundary()[patchi].edgeFaces();

        if (boundary()[patchi].coupled())
        {
            vectorField nbrCf(Cf.boundaryField()[patchi].patchNeighbourField());
            vectorField nbrNf(Nf.boundaryField()[patchi].patchNeighbourField());

            forAll(edgeFaces, bndEdgei)
            {
                const label ownFacei = edgeFaces[bndEdgei];
                const label meshEdgei(boundary()[patchi].start() + bndEdgei);

                auto& tensors = edgeTransformTensors[meshEdgei];

                vector edgeCtr = Ce.boundaryField()[patchi][bndEdgei];

                if (skew())
                {
                    edgeCtr -= skewCorrectionVectors()
                        .boundaryField()[patchi][bndEdgei];
                }

                // Edge transformation tensor
                tensors[0] =
                    rotation_e3e1
                    (
                        Ne.boundaryField()[patchi][bndEdgei],
                        (edgeCtr - Cf.internalField()[ownFacei])
                    );

                // Owner transformation tensor
                tensors[1] =
                    rotation_e3e1
                    (
                        Nf.internalField()[ownFacei],
                        (edgeCtr - Cf.internalField()[ownFacei])
                    );

                // Neighbour transformation tensor
                tensors[2] =
                    rotation_e3e1
                    (
                        nbrNf[bndEdgei],
                        (nbrCf[bndEdgei] - edgeCtr)
                    );
            }
        }
        else
        {
            forAll(edgeFaces, bndEdgei)
            {
                const label ownFacei = edgeFaces[bndEdgei];
                const label meshEdgei(boundary()[patchi].start() + bndEdgei);

                auto& tensors = edgeTransformTensors[meshEdgei];

                vector edgeCtr = Ce.boundaryField()[patchi][bndEdgei];

                if (skew())
                {
                    edgeCtr -= skewCorrectionVectors()
                        .boundaryField()[patchi][bndEdgei];
                }

                // Edge transformation tensor
                tensors[0] =
                    rotation_e3e1
                    (
                        Ne.boundaryField()[patchi][bndEdgei],
                        (edgeCtr - Cf.internalField()[ownFacei])
                    );

                // Owner transformation tensor
                tensors[1] =
                    rotation_e3e1
                    (
                        Nf.internalField()[ownFacei],
                        (edgeCtr - Cf.internalField()[ownFacei])
                    );

                // Neighbour transformation tensor
                tensors[2] = tensor::I;
            }
        }
    }
}


Foam::labelList Foam::faMesh::internalPoints() const
{
    DebugInFunction
        << "Calculating internal points" << endl;

    bitSet markPoints(markupBoundaryPoints(this->patch()));
    markPoints.flip();

    return markPoints.sortedToc();
}


Foam::labelList Foam::faMesh::boundaryPoints() const
{
    DebugInFunction
        << "Calculating boundary points" << endl;

    bitSet markPoints(markupBoundaryPoints(this->patch()));

    return markPoints.sortedToc();
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~
// Point normal calculations
// ~~~~~~~~~~~~~~~~~~~~~~~~~

// Revised method (general)
// ------------------------
//
// - For each patch face obtain a centre point (mathematical avg)
//   and use that to define the face dual as a pair of triangles:
//     - tri1: base-point / mid-side of right edge / face centre
//     - tri2: base-point / face centre / mid-side of left edge
//
// - Walk all face points, inserting directly into the corresponding
//   locations. No distinction between internal or boundary points (yet).
//
// Revised method (boundary correction)
// ------------------------------------
//
// - correct wedge directly, use processor patch information to exchange
//   the current summed values. [No change from original].
//
// - explicit correction of other boundaries.
//   Use the new boundary halo information for the face normals.
//   Calculate the equivalent face-point normals locally and apply
//   correction as before.

void Foam::faMesh::calcPointAreaNormals(vectorField& result) const
{
    DebugInFunction
        << "Calculating pointAreaNormals : face dual method" << endl;

    result.resize_nocopy(nPoints());
    result = Zero;

    const pointField& points = patch().localPoints();
    const faceList& faces = patch().localFaces();


    // Loop over all faces

    for (const face& f : faces)
    {
        const label nVerts(f.size());

        point centrePoint(Zero);
        for (label i = 0; i < nVerts; ++i)
        {
            centrePoint += points[f[i]];
        }
        centrePoint /= nVerts;

        for (label i = 0; i < nVerts; ++i)
        {
            const label pt0 = f.thisLabel(i);    // base

            result[pt0] +=
                areaInvDistSqrWeightedNormalDualEdge
                (
                    points[pt0],             // base
                    points[f.nextLabel(i)],  // right
                    points[f.prevLabel(i)],  // left
                    centrePoint
                );
        }
    }


    // Handle the boundary edges

    bitSet nbrBoundaryAdjust(boundary().size(), true);

    forAll(boundary(), patchi)
    {
        const faPatch& fap = boundary()[patchi];

        if (isA<wedgeFaPatch>(fap))
        {
            // Correct wedge points

            const auto& wedgePatch = refCast<const wedgeFaPatch>(fap);

            const labelList& patchPoints = wedgePatch.pointLabels();

            const vector N
            (
                transform
                (
                    wedgePatch.edgeT(),
                    wedgePatch.centreNormal()
                ).normalise()
            );

            for (const label pti : patchPoints)
            {
                result[pti].removeCollinear(N);
            }

            // Axis point correction
            if (wedgePatch.axisPoint() > -1)
            {
                result[wedgePatch.axisPoint()] =
                    wedgePatch.axis()
                   *(
                       wedgePatch.axis()
                      &result[wedgePatch.axisPoint()]
                    );
            }

            // Handled
            nbrBoundaryAdjust.unset(patchi);
        }
        else if (Pstream::parRun() && isA<processorFaPatch>(fap))
        {
            // Correct processor patch points

            const auto& procPatch = refCast<const processorFaPatch>(fap);

            const labelList& patchPoints = procPatch.pointLabels();
            const labelList& nbrPatchPoints = procPatch.neighbPoints();

            const labelList& nonGlobalPatchPoints =
                procPatch.nonGlobalPatchPoints();

            // Send my values

            vectorField patchPointNormals
            (
                UIndirectList<vector>(result, patchPoints)
            );

            {
                OPstream::write
                (
                    Pstream::commsTypes::blocking,
                    procPatch.neighbProcNo(),
                    patchPointNormals.cdata_bytes(),
                    patchPointNormals.size_bytes()
                );
            }

            // Receive neighbour values
            patchPointNormals.resize(nbrPatchPoints.size());

            {
                IPstream::read
                (
                    Pstream::commsTypes::blocking,
                    procPatch.neighbProcNo(),
                    patchPointNormals.data_bytes(),
                    patchPointNormals.size_bytes()
                );
            }

            for (const label pti : nonGlobalPatchPoints)
            {
                result[patchPoints[pti]] +=
                    patchPointNormals[nbrPatchPoints[pti]];
            }

            // Handled
            nbrBoundaryAdjust.unset(patchi);
        }
        else if (fap.coupled())
        {
            // Coupled - no further action for neighbour side
            nbrBoundaryAdjust.unset(patchi);
        }
        // TBD:
        /// else if (isA<emptyFaPatch>(fap))
        /// {
        ///     // Ignore this boundary
        ///     nbrBoundaryAdjust.unset(patchi);
        /// }
        else if (!correctPatchPointNormals(patchi))
        {
            // No corrections
            nbrBoundaryAdjust.unset(patchi);
        }
    }


    // Correct global points
    if (globalData().nGlobalPoints())
    {
        const labelList& spLabels = globalData().sharedPointLabels();
        const labelList& addr = globalData().sharedPointAddr();

        vectorField spNormals
        (
            UIndirectList<vector>(result, spLabels)
        );

        vectorField gpNormals
        (
            globalData().nGlobalPoints(),
            Zero
        );

        forAll(addr, i)
        {
            gpNormals[addr[i]] += spNormals[i];
        }

        Pstream::combineAllGather(gpNormals, plusEqOp<vectorField>());

        // Extract local data
        forAll(addr, i)
        {
            spNormals[i] = gpNormals[addr[i]];
        }

        forAll(spNormals, pointI)
        {
            result[spLabels[pointI]] = spNormals[pointI];
        }
    }


    if (returnReduce(nbrBoundaryAdjust.any(), orOp<bool>()))
    {
        if (debug)
        {
            PoutInFunction
                << "Apply " << nbrBoundaryAdjust.count()
                << " boundary neighbour corrections" << nl;
        }

        // Apply boundary points correction

        // Collect face normals as point normals

        const auto& haloNormals = this->haloFaceNormals();

        Map<vector> fpNormals(4*nBoundaryEdges());

        for (const label patchi : nbrBoundaryAdjust)
        {
            const faPatch& fap = boundary()[patchi];
            const labelList& edgeLabels = fap.edgeLabels();

            if (fap.ngbPolyPatchIndex() < 0)
            {
                FatalErrorInFunction
                    << "Neighbour polyPatch index is not defined "
                    << "for faPatch " << fap.name()
                    << abort(FatalError);
            }

            for (const label edgei : edgeLabels)
            {
                const edge& e = patch().edges()[edgei];

                // Halo face unitNormal at boundary edge (starts as 0)
                const vector& fnorm = haloNormals[edgei - nInternalEdges_];

                fpNormals(e.first()) += fnorm;
                fpNormals(e.second()) += fnorm;
            }
        }

        // Apply the correction

        // Note from Zeljko Tukovic:
        //
        // This posibility is used for free-surface tracking
        // calculations to enforce 90 deg contact angle between
        // finite-area mesh and neighbouring polyPatch. It is very
        // important for curvature calculation to have correct normal
        // at contact line points.

        forAllConstIters(fpNormals, iter)
        {
            const label pointi = iter.key();
            vector fpnorm = normalised(iter.val());

            result[pointi].removeCollinear(fpnorm);
        }
    }

    result.normalise();
}


void Foam::faMesh::calcPointAreaNormalsByQuadricsFit(vectorField& result) const
{
    const labelList intPoints(internalPoints());
    const labelList bndPoints(boundaryPoints());

    const pointField& points = patch().localPoints();
    const faceList& faces = patch().localFaces();
    const labelListList& pointFaces = patch().pointFaces();

    forAll(intPoints, pointI)
    {
        label curPoint = intPoints[pointI];

        const labelHashSet faceSet(pointFaces[curPoint]);
        const labelList curFaces(faceSet.toc());

        labelHashSet pointSet;

        for (const label facei : curFaces)
        {
            const labelList& facePoints = faces[facei];
            pointSet.insert(facePoints);
        }
        pointSet.erase(curPoint);
        labelList curPoints(pointSet.toc());

        if (pointSet.size() < 5)
        {
            DebugInfo
                << "WARNING: Extending point set for fitting." << endl;

            labelHashSet faceSet(pointFaces[curPoint]);
            labelList curFaces(faceSet.toc());
            for (const label facei : curFaces)
            {
                const labelList& curFaceFaces = patch().faceFaces()[facei];
                faceSet.insert(curFaceFaces);
            }
            curFaces = faceSet.toc();

            pointSet.clear();

            for (const label facei : curFaces)
            {
                const labelList& facePoints = faces[facei];
                pointSet.insert(facePoints);
            }
            pointSet.erase(curPoint);
            curPoints = pointSet.toc();
        }

        pointField allPoints(curPoints.size());
        scalarField W(curPoints.size(), 1.0);
        for (label i=0; i<curPoints.size(); ++i)
        {
            allPoints[i] = points[curPoints[i]];
            W[i] = 1.0/magSqr(allPoints[i] - points[curPoint]);
        }

        // Transform points
        coordSystem::cartesian cs
        (
            points[curPoint],   // origin
            result[curPoint],   // axis [e3] (normalized by constructor)
            allPoints[0] - points[curPoint] // direction [e1]
        );

        for (point& p : allPoints)
        {
            p = cs.localPosition(p);
        }

        scalarRectangularMatrix M
        (
            allPoints.size(),
            5,
            0.0
        );

        for (label i = 0; i < allPoints.size(); ++i)
        {
            M[i][0] = sqr(allPoints[i].x());
            M[i][1] = sqr(allPoints[i].y());
            M[i][2] = allPoints[i].x()*allPoints[i].y();
            M[i][3] = allPoints[i].x();
            M[i][4] = allPoints[i].y();
        }

        scalarSquareMatrix MtM(5, Zero);

        for (label i = 0; i < MtM.n(); ++i)
        {
            for (label j = 0; j < MtM.m(); ++j)
            {
                for (label k = 0; k < M.n(); ++k)
                {
                    MtM[i][j] += M[k][i]*M[k][j]*W[k];
                }
            }
        }

        scalarField MtR(5, Zero);

        for (label i=0; i<MtR.size(); ++i)
        {
            for (label j=0; j<M.n(); ++j)
            {
                MtR[i] += M[j][i]*allPoints[j].z()*W[j];
            }
        }

        LUsolve(MtM, MtR);

        vector curNormal = vector(MtR[3], MtR[4], -1);

        curNormal = cs.globalVector(curNormal);

        curNormal *= sign(curNormal&result[curPoint]);

        result[curPoint] = curNormal;
    }


    forAll(boundary(), patchI)
    {
        const faPatch& fap = boundary()[patchI];

        if (Pstream::parRun() && isA<processorFaPatch>(fap))
        {
            const processorFaPatch& procPatch =
                refCast<const processorFaPatch>(boundary()[patchI]);

            const labelList& patchPointLabels = procPatch.pointLabels();

            labelList toNgbProcLsPointStarts(patchPointLabels.size(), Zero);
            vectorField toNgbProcLsPoints
            (
                10*patchPointLabels.size(),
                Zero
            );
            label nPoints = 0;

            for (label pointI=0; pointI<patchPointLabels.size(); ++pointI)
            {
                label curPoint = patchPointLabels[pointI];

                toNgbProcLsPointStarts[pointI] = nPoints;

                const labelHashSet faceSet(pointFaces[curPoint]);
                const labelList curFaces(faceSet.toc());

                labelHashSet pointSet;

                for (const label facei : curFaces)
                {
                    const labelList& facePoints = faces[facei];
                    pointSet.insert(facePoints);
                }
                pointSet.erase(curPoint);
                labelList curPoints = pointSet.toc();

                for (label i=0; i<curPoints.size(); ++i)
                {
                    toNgbProcLsPoints[nPoints++] = points[curPoints[i]];
                }
            }

            toNgbProcLsPoints.setSize(nPoints);

            {
                OPstream toNeighbProc
                (
                    Pstream::commsTypes::blocking,
                    procPatch.neighbProcNo(),
                    toNgbProcLsPoints.size_bytes()
                  + toNgbProcLsPointStarts.size_bytes()
                  + 10*sizeof(label)
                );

                toNeighbProc
                    << toNgbProcLsPoints
                    << toNgbProcLsPointStarts;
            }
        }
    }

    for (const faPatch& fap : boundary())
    {
        if (Pstream::parRun() && isA<processorFaPatch>(fap))
        {
            const processorFaPatch& procPatch =
                refCast<const processorFaPatch>(fap);

            const labelList& patchPointLabels = procPatch.pointLabels();

            labelList fromNgbProcLsPointStarts(patchPointLabels.size(), Zero);
            vectorField fromNgbProcLsPoints;

            {
                IPstream fromNeighbProc
                (
                    Pstream::commsTypes::blocking,
                    procPatch.neighbProcNo(),
                    10*patchPointLabels.size()*sizeof(vector)
                  + fromNgbProcLsPointStarts.size_bytes()
                  + 10*sizeof(label)
                );

                fromNeighbProc
                    >> fromNgbProcLsPoints
                    >> fromNgbProcLsPointStarts;
            }

            const labelList& nonGlobalPatchPoints =
                procPatch.nonGlobalPatchPoints();

            forAll(nonGlobalPatchPoints, pointI)
            {
                label curPoint =
                    patchPointLabels[nonGlobalPatchPoints[pointI]];
                label curNgbPoint =
                    procPatch.neighbPoints()[nonGlobalPatchPoints[pointI]];

                labelHashSet faceSet;
                faceSet.insert(pointFaces[curPoint]);

                labelList curFaces = faceSet.toc();

                labelHashSet pointSet;

                for (const label facei : curFaces)
                {
                    const labelList& facePoints = faces[facei];
                    pointSet.insert(facePoints);
                }
                pointSet.erase(curPoint);
                labelList curPoints = pointSet.toc();

                label nAllPoints = curPoints.size();

                if (curNgbPoint == fromNgbProcLsPointStarts.size() - 1)
                {
                    nAllPoints +=
                        fromNgbProcLsPoints.size()
                      - fromNgbProcLsPointStarts[curNgbPoint];
                }
                else
                {
                    nAllPoints +=
                        fromNgbProcLsPointStarts[curNgbPoint + 1]
                      - fromNgbProcLsPointStarts[curNgbPoint];
                }

                vectorField allPointsExt(nAllPoints);
                label counter = 0;
                for (label i=0; i<curPoints.size(); ++i)
                {
                    allPointsExt[counter++] = points[curPoints[i]];
                }

                if (curNgbPoint == fromNgbProcLsPointStarts.size() - 1)
                {
                    for
                    (
                        label i=fromNgbProcLsPointStarts[curNgbPoint];
                        i<fromNgbProcLsPoints.size();
                        ++i
                    )
                    {
                        allPointsExt[counter++] = fromNgbProcLsPoints[i];
                    }
                }
                else
                {
                    for
                    (
                        label i=fromNgbProcLsPointStarts[curNgbPoint];
                        i<fromNgbProcLsPointStarts[curNgbPoint+1];
                        ++i
                    )
                    {
                        allPointsExt[counter++] = fromNgbProcLsPoints[i];
                    }
                }

                // Remove duplicate points
                vectorField allPoints(nAllPoints, Zero);
                boundBox bb(allPointsExt, false);
                scalar tol = 0.001*mag(bb.max() - bb.min());

                nAllPoints = 0;
                forAll(allPointsExt, pI)
                {
                    bool duplicate = false;
                    for (label i=0; i<nAllPoints; ++i)
                    {
                        if
                        (
                            mag
                            (
                                allPoints[i]
                              - allPointsExt[pI]
                            )
                          < tol
                        )
                        {
                            duplicate = true;
                            break;
                        }
                    }

                    if (!duplicate)
                    {
                        allPoints[nAllPoints++] = allPointsExt[pI];
                    }
                }

                allPoints.setSize(nAllPoints);

                if (nAllPoints < 5)
                {
                    FatalErrorInFunction
                        << "There are no enough points for quadrics "
                        << "fitting for a point at processor patch"
                        << abort(FatalError);
                }

                // Transform points
                const vector& origin = points[curPoint];
                const vector axis = normalised(result[curPoint]);
                vector dir(allPoints[0] - points[curPoint]);
                dir.removeCollinear(axis);
                dir.normalise();

                const coordinateSystem cs("cs", origin, axis, dir);

                scalarField W(allPoints.size(), 1.0);

                forAll(allPoints, pI)
                {
                    W[pI] = 1.0/magSqr(allPoints[pI] - points[curPoint]);

                    allPoints[pI] = cs.localPosition(allPoints[pI]);
                }

                scalarRectangularMatrix M
                (
                    allPoints.size(),
                    5,
                    0.0
                );

                for (label i=0; i<allPoints.size(); ++i)
                {
                    M[i][0] = sqr(allPoints[i].x());
                    M[i][1] = sqr(allPoints[i].y());
                    M[i][2] = allPoints[i].x()*allPoints[i].y();
                    M[i][3] = allPoints[i].x();
                    M[i][4] = allPoints[i].y();
                }

                scalarSquareMatrix MtM(5, Zero);

                for (label i = 0; i < MtM.n(); ++i)
                {
                    for (label j = 0; j < MtM.m(); ++j)
                    {
                        for (label k = 0; k < M.n(); ++k)
                        {
                            MtM[i][j] += M[k][i]*M[k][j]*W[k];
                        }
                    }
                }

                scalarField MtR(5, Zero);

                for (label i = 0; i < MtR.size(); ++i)
                {
                    for (label j = 0; j < M.n(); ++j)
                    {
                        MtR[i] += M[j][i]*allPoints[j].z()*W[j];
                    }
                }

                LUsolve(MtM, MtR);

                vector curNormal = vector(MtR[3], MtR[4], -1);

                curNormal = cs.globalVector(curNormal);

                curNormal *= sign(curNormal&result[curPoint]);

                result[curPoint] = curNormal;
            }
        }
    }

    // Correct global points
    if (globalData().nGlobalPoints() > 0)
    {
        const labelList& spLabels = globalData().sharedPointLabels();

        const labelList& addr = globalData().sharedPointAddr();

        for (label k=0; k<globalData().nGlobalPoints(); ++k)
        {
            List<List<vector>> procLsPoints(Pstream::nProcs());

            const label curSharedPointIndex = addr.find(k);

            scalar tol = 0.0;

            if (curSharedPointIndex != -1)
            {
                label curPoint = spLabels[curSharedPointIndex];

                const labelHashSet faceSet(pointFaces[curPoint]);
                const labelList curFaces(faceSet.toc());

                labelHashSet pointSet;
                for (const label facei : curFaces)
                {
                    const labelList& facePoints = faces[facei];
                    pointSet.insert(facePoints);
                }
                pointSet.erase(curPoint);
                labelList curPoints = pointSet.toc();

                vectorField locPoints(points, curPoints);

                procLsPoints[Pstream::myProcNo()] = locPoints;

                const boundBox bb(locPoints, false);
                tol = 0.001*mag(bb.max() - bb.min());
            }

            Pstream::allGatherList(procLsPoints);

            if (curSharedPointIndex != -1)
            {
                label curPoint = spLabels[curSharedPointIndex];

                label nAllPoints = 0;
                forAll(procLsPoints, procI)
                {
                    nAllPoints += procLsPoints[procI].size();
                }

                vectorField allPoints(nAllPoints, Zero);

                nAllPoints = 0;
                forAll(procLsPoints, procI)
                {
                    forAll(procLsPoints[procI], pointI)
                    {
                        bool duplicate = false;
                        for (label i=0; i<nAllPoints; ++i)
                        {
                            if
                            (
                                mag
                                (
                                    allPoints[i]
                                  - procLsPoints[procI][pointI]
                                )
                              < tol
                            )
                            {
                                duplicate = true;
                                break;
                            }
                        }

                        if (!duplicate)
                        {
                            allPoints[nAllPoints++] =
                                procLsPoints[procI][pointI];
                        }
                    }
                }

                allPoints.setSize(nAllPoints);

                if (nAllPoints < 5)
                {
                    FatalErrorInFunction
                        << "There are no enough points for quadratic "
                        << "fitting for a global processor point "
                        << abort(FatalError);
                }

                // Transform points
                const vector& origin = points[curPoint];
                const vector axis = normalised(result[curPoint]);
                vector dir(allPoints[0] - points[curPoint]);
                dir.removeCollinear(axis);
                dir.normalise();

                coordinateSystem cs("cs", origin, axis, dir);

                scalarField W(allPoints.size(), 1.0);

                forAll(allPoints, pointI)
                {
                    W[pointI]=
                        1.0/magSqr(allPoints[pointI] - points[curPoint]);

                    allPoints[pointI] = cs.localPosition(allPoints[pointI]);
                }

                scalarRectangularMatrix M
                (
                    allPoints.size(),
                    5,
                    0.0
                );

                for (label i=0; i<allPoints.size(); ++i)
                {
                    M[i][0] = sqr(allPoints[i].x());
                    M[i][1] = sqr(allPoints[i].y());
                    M[i][2] = allPoints[i].x()*allPoints[i].y();
                    M[i][3] = allPoints[i].x();
                    M[i][4] = allPoints[i].y();
                }

                scalarSquareMatrix MtM(5, Zero);
                for (label i = 0; i < MtM.n(); ++i)
                {
                    for (label j = 0; j < MtM.m(); ++j)
                    {
                        for (label k = 0; k < M.n(); ++k)
                        {
                            MtM[i][j] += M[k][i]*M[k][j]*W[k];
                        }
                    }
                }

                scalarField MtR(5, Zero);
                for (label i = 0; i < MtR.size(); ++i)
                {
                    for (label j = 0; j < M.n(); ++j)
                    {
                        MtR[i] += M[j][i]*allPoints[j].z()*W[j];
                    }
                }

                LUsolve(MtM, MtR);

                vector curNormal = vector(MtR[3], MtR[4], -1);

                curNormal = cs.globalVector(curNormal);

                curNormal *= sign(curNormal&result[curPoint]);

                result[curPoint] = curNormal;
            }
        }
    }

    for (vector& n : result)
    {
        n.normalise();
    }
}


Foam::tmp<Foam::edgeScalarField> Foam::faMesh::edgeLengthCorrection() const
{
    DebugInFunction
        << "Calculating edge length correction" << endl;

    auto tcorrection = tmp<edgeScalarField>::New
    (
        IOobject
        (
            "edgeLengthCorrection",
            mesh().pointsInstance(),
            meshSubDir,
            mesh()
        ),
        *this,
        dimless
    );
    auto& correction = tcorrection.ref();

    const vectorField& pointNormals = pointAreaNormals();

    forAll(correction.internalField(), edgeI)
    {
        scalar sinAlpha = mag
        (
            pointNormals[edges()[edgeI].start()]^
            pointNormals[edges()[edgeI].end()]
        );

        scalar alpha = asin(sinAlpha);

        correction.ref()[edgeI] = cos(0.5*alpha);
    }


    forAll(boundary(), patchI)
    {
        const edgeList::subList patchEdges
        (
             boundary()[patchI].patchSlice(edges())
        );

        forAll(patchEdges, edgeI)
        {
            scalar sinAlpha =
                mag
                (
                    pointNormals[patchEdges[edgeI].start()]
                  ^ pointNormals[patchEdges[edgeI].end()]
                );

            scalar alpha = asin(sinAlpha);

            correction.boundaryFieldRef()[patchI][edgeI] = cos(0.5*alpha);
        }
    }

    return tcorrection;
}


// ************************************************************************* //
