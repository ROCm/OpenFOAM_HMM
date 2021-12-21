/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
#include "primitiveFacePatch.H"
#include "fac.H"
#include "processorFaPatch.H"
#include "wedgeFaPatch.H"
#include "PstreamCombineReduceOps.H"
#include "cartesianCS.H"
#include "scalarMatrices.H"
#include "processorFaPatchFields.H"
#include "emptyFaPatchFields.H"

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


    const pointField& pPoints = points();
    const edgeList& pEdges = edges();

    const edgeVectorField& eCentres = edgeCentres();
    const areaVectorField& fCentres = areaCentres();

    const edgeVectorField& edgeNormals = edgeAreaNormals();

    vectorField& leInternal = Le.ref();
    const vectorField& edgeNormalsInternal = edgeNormals.internalField();
    const vectorField& fCentresInternal = fCentres.internalField();
    const vectorField& eCentresInternal = eCentres.internalField();
    const scalarField& magLeInternal = magLe().internalField();

    forAll(leInternal, edgeI)
    {
        leInternal[edgeI] =
            pEdges[edgeI].vec(pPoints) ^ edgeNormalsInternal[edgeI];

        leInternal[edgeI] *=
          - sign
            (
                leInternal[edgeI] &
                (
                    fCentresInternal[owner()[edgeI]]
                  - eCentresInternal[edgeI]
                )
            );

        leInternal[edgeI] *=
            magLeInternal[edgeI]/mag(leInternal[edgeI]);
    }

    forAll(boundary(), patchI)
    {
        const labelUList& bndEdgeFaces = boundary()[patchI].edgeFaces();

        const edgeList::subList bndEdges =
            boundary()[patchI].patchSlice(pEdges);

        const vectorField& bndEdgeNormals =
            edgeNormals.boundaryField()[patchI];

        vectorField& patchLe = Le.boundaryFieldRef()[patchI];
        const vectorField& patchECentres = eCentres.boundaryField()[patchI];

        forAll(patchLe, edgeI)
        {
            patchLe[edgeI] =
                bndEdges[edgeI].vec(pPoints) ^ bndEdgeNormals[edgeI];

            patchLe[edgeI] *=
              - sign
                (
                    patchLe[edgeI]&
                    (
                        fCentresInternal[bndEdgeFaces[edgeI]]
                      - patchECentres[edgeI]
                    )
                );

            patchLe[edgeI] *=
                magLe().boundaryField()[patchI][edgeI]
                /mag(patchLe[edgeI]);
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

    label edgei = 0;
    for (const edge& e : patch().internalEdges())
    {
        magLe.ref()[edgei] = e.mag(localPoints);
        ++edgei;
    }


    forAll(boundary(), patchI)
    {
        const edgeList::subList patchEdges =
            boundary()[patchI].patchSlice(edges());

        forAll(patchEdges, edgeI)
        {
            magLe.boundaryFieldRef()[patchI][edgeI] =
                patchEdges[edgeI].mag(localPoints);
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

    forAll(localFaces, faceI)
    {
        centres.ref()[faceI] = localFaces[faceI].centre(localPoints);
    }

    forAll(boundary(), patchI)
    {
        const edgeList::subList patchEdges =
            boundary()[patchI].patchSlice(edges());

        forAll(patchEdges, edgeI)
        {
            centres.boundaryFieldRef()[patchI][edgeI] =
                patchEdges[edgeI].centre(localPoints);
        }
    }

    forAll(centres.boundaryField(), patchI)
    {
        //HJ: this is wrong!  5/Aug/2011
        if
        (
            isA<processorFaPatchVectorField>
            (
                centres.boundaryField()[patchI]
            )
        )
        {
            centres.boundaryFieldRef()[patchI].initEvaluate();
            centres.boundaryFieldRef()[patchI].evaluate();
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

    edgeVectorField& edgeCentres = *edgeCentresPtr_;

    const pointField& localPoints = points();


    label edgei = 0;
    for (const edge& e : patch().internalEdges())
    {
        edgeCentres.ref()[edgei] = e.centre(localPoints);
        ++edgei;
    }


    forAll(boundary(), patchI)
    {
        const edgeList::subList patchEdges =
            boundary()[patchI].patchSlice(edges());

        forAll(patchEdges, edgeI)
        {
            edgeCentres.boundaryFieldRef()[patchI][edgeI] =
                patchEdges[edgeI].centre(localPoints);
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
    DimensionedField<scalar, areaMesh>& S = *SPtr_;

    const pointField& localPoints = points();
    const faceList& localFaces = faces();

    forAll(S, faceI)
    {
        S[faceI] = localFaces[faceI].mag(localPoints);
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

    areaVectorField& faceAreaNormals = *faceAreaNormalsPtr_;

    const pointField& localPoints = points();
    const faceList& localFaces = faces();

    vectorField& nInternal = faceAreaNormals.ref();
    forAll(localFaces, faceI)
    {
        nInternal[faceI] = localFaces[faceI].unitNormal(localPoints);
    }

    forAll(boundary(), patchI)
    {
        faceAreaNormals.boundaryFieldRef()[patchI] =
            edgeAreaNormals().boundaryField()[patchI];
    }

    forAll(faceAreaNormals.boundaryField(), patchI)
    {
        if
        (
            isA<processorFaPatchVectorField>
            (
                faceAreaNormals.boundaryField()[patchI]
            )
        )
        {
            faceAreaNormals.boundaryFieldRef()[patchI].initEvaluate();
            faceAreaNormals.boundaryFieldRef()[patchI].evaluate();
        }
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


    // Point area normals
    const vectorField& pointNormals = pointAreaNormals();

    forAll(edgeAreaNormals.internalField(), edgei)
    {
        const edge& e = edges()[edgei];
        const vector edgeVec = e.unitVec(points());

        vector& n = edgeAreaNormals.ref()[edgei];

        n = (pointNormals[e.first()] + pointNormals[e.second()]);

        n -= edgeVec*(edgeVec & n);

        n.normalise();
    }

    forAll(boundary(), patchi)
    {
        const edgeList::subList patchEdges =
            boundary()[patchi].patchSlice(edges());

        vectorField& edgeNorms = edgeAreaNormals.boundaryFieldRef()[patchi];

        forAll(patchEdges, edgei)
        {
            const edge& e = patchEdges[edgei];
            const vector edgeVec = e.unitVec(points());

            vector& n = edgeNorms[edgei];

            n = (pointNormals[e.first()] + pointNormals[e.second()]);

            n -= edgeVec*(edgeVec & n);

            n.normalise();
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

    edgeTransformTensorsPtr_ = new FieldField<Field, tensor>(nEdges());
    FieldField<Field, tensor>& edgeTransformTensors =
        *edgeTransformTensorsPtr_;

    const areaVectorField& Nf = faceAreaNormals();
    const areaVectorField& Cf = areaCentres();

    const edgeVectorField& Ne = edgeAreaNormals();
    const edgeVectorField& Ce = edgeCentres();

    // Internal edges transformation tensors
    for (label edgeI=0; edgeI<nInternalEdges(); ++edgeI)
    {
        edgeTransformTensors.set(edgeI, new Field<tensor>(3, I));

        vector E = Ce.internalField()[edgeI];

        if (skew())
        {
            E -= skewCorrectionVectors().internalField()[edgeI];
        }

        // Edge transformation tensor
        vector il = E - Cf.internalField()[owner()[edgeI]];

        il -= Ne.internalField()[edgeI]
            *(Ne.internalField()[edgeI]&il);

        il /= mag(il);

        vector kl = Ne.internalField()[edgeI];
        vector jl = kl^il;

        edgeTransformTensors[edgeI][0] =
            tensor
            (
                il.x(), il.y(), il.z(),
                jl.x(), jl.y(), jl.z(),
                kl.x(), kl.y(), kl.z()
            );

        // Owner transformation tensor
        il = E - Cf.internalField()[owner()[edgeI]];

        il -= Nf.internalField()[owner()[edgeI]]
            *(Nf.internalField()[owner()[edgeI]]&il);

        il /= mag(il);

        kl = Nf.internalField()[owner()[edgeI]];
        jl = kl^il;

        edgeTransformTensors[edgeI][1] =
            tensor
            (
                il.x(), il.y(), il.z(),
                jl.x(), jl.y(), jl.z(),
                kl.x(), kl.y(), kl.z()
            );

        // Neighbour transformation tensor
        il = Cf.internalField()[neighbour()[edgeI]] - E;

        il -= Nf.internalField()[neighbour()[edgeI]]
            *(Nf.internalField()[neighbour()[edgeI]]&il);

        il /= mag(il);

        kl = Nf.internalField()[neighbour()[edgeI]];
        jl = kl^il;

        edgeTransformTensors[edgeI][2] =
            tensor
            (
                il.x(), il.y(), il.z(),
                jl.x(), jl.y(), jl.z(),
                kl.x(), kl.y(), kl.z()
            );
    }

    // Boundary edges transformation tensors
    forAll(boundary(), patchI)
    {
        if (boundary()[patchI].coupled())
        {
            const labelUList& edgeFaces =
                boundary()[patchI].edgeFaces();

            vectorField ngbCf(Cf.boundaryField()[patchI].patchNeighbourField());

            vectorField ngbNf(Nf.boundaryField()[patchI].patchNeighbourField());

            forAll(edgeFaces, edgeI)
            {
                edgeTransformTensors.set
                (
                    boundary()[patchI].start() + edgeI,
                    new Field<tensor>(3, I)
                );

                vector E = Ce.boundaryField()[patchI][edgeI];

                if (skew())
                {
                    E -= skewCorrectionVectors()
                        .boundaryField()[patchI][edgeI];
                }

                // Edge transformation tensor
                vector il = E - Cf.internalField()[edgeFaces[edgeI]];

                il -= Ne.boundaryField()[patchI][edgeI]
                   *(Ne.boundaryField()[patchI][edgeI]&il);

                il /= mag(il);

                vector kl = Ne.boundaryField()[patchI][edgeI];
                vector jl = kl^il;

                edgeTransformTensors[boundary()[patchI].start() + edgeI][0] =
                    tensor(il, jl, kl);

                // Owner transformation tensor
                il = E - Cf.internalField()[edgeFaces[edgeI]];

                il -= Nf.internalField()[edgeFaces[edgeI]]
                   *(Nf.internalField()[edgeFaces[edgeI]]&il);

                il /= mag(il);

                kl = Nf.internalField()[edgeFaces[edgeI]];
                jl = kl^il;

                edgeTransformTensors[boundary()[patchI].start() + edgeI][1] =
                    tensor(il, jl, kl);

                // Neighbour transformation tensor
                il = ngbCf[edgeI] - E;

                il -= ngbNf[edgeI]*(ngbNf[edgeI]&il);

                il /= mag(il);

                kl = ngbNf[edgeI];

                jl = kl^il;

                edgeTransformTensors[boundary()[patchI].start() + edgeI][2] =
                    tensor(il, jl, kl);
            }
        }
        else
        {
            const labelUList& edgeFaces = boundary()[patchI].edgeFaces();

            forAll(edgeFaces, edgeI)
            {
                edgeTransformTensors.set
                (
                    boundary()[patchI].start() + edgeI,
                    new Field<tensor>(3, I)
                );

                vector E = Ce.boundaryField()[patchI][edgeI];

                if (skew())
                {
                    E -= skewCorrectionVectors()
                        .boundaryField()[patchI][edgeI];
                }

                // Edge transformation tensor
                vector il = E - Cf.internalField()[edgeFaces[edgeI]];

                il -= Ne.boundaryField()[patchI][edgeI]
                   *(Ne.boundaryField()[patchI][edgeI]&il);

                il /= mag(il);

                vector kl = Ne.boundaryField()[patchI][edgeI];
                vector jl = kl^il;

                edgeTransformTensors[boundary()[patchI].start() + edgeI][0] =
                    tensor(il, jl, kl);

                // Owner transformation tensor
                il = E - Cf.internalField()[edgeFaces[edgeI]];

                il -= Nf.internalField()[edgeFaces[edgeI]]
                   *(Nf.internalField()[edgeFaces[edgeI]]&il);

                il /= mag(il);

                kl = Nf.internalField()[edgeFaces[edgeI]];
                jl = kl^il;

                edgeTransformTensors[boundary()[patchI].start() + edgeI][1] =
                    tensor(il, jl, kl);
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

// Original method (general)
// -------------------------
// - For each point, obtain the list of connected patch faces
//   (from point-to-face addressing).
//
// - Create a primitive patch for those faces and use that to obtain the
//   outer edge loop(s). This is effectively an agglomeration of the patch
//   faces connected to a point.
//
// - Perform a pair-wise walk of the edge loop to obtain triangles from
//   the originating point outwards (fan-like triangulation).
//   Calculate an area-weighted value for each triangle.
//
//   NOTE: not sure why internal and boundary point agglomeration was
//   handled separately.
//
// Problems:
// - possibly susceptible to edge-loop errors (issue #2233) that cause
//   the agglomeration logic to include the current point twice?
// - use of outer edge loop makes it more sensitive to face warpage.
// - relatively expensive with point-face connectivity,
//   creation/destruction of a primitive-patch around each point.
//
// Original method (boundary correction)
// -------------------------------------
//
// - correct wedge directly, use processor patch information to exchange
//   the current summed values
//
// - explicit correction of other boundaries.
//   Polls the patch for the ngbPolyPatchPointNormals(), which internally
//   calls ngbPolyPatchFaces and can return -1 for unmatched edges.
//   This occurs when the outside perimeter of the faPatch aligns with
//   a polyMesh processor. The neighbour face is off-processor and cannot
//   be found. Accessing the mesh face at -1 == SEGFAULT.

void Foam::faMesh::calcPointAreaNormals_orig(vectorField& result) const
{
    DebugInFunction
        << "Calculating pointAreaNormals : original method" << endl;

    result.resize(nPoints());
    result = Zero;

    labelList intPoints(internalPoints());
    labelList bndPoints(boundaryPoints());

    const pointField& points = patch().localPoints();
    const faceList& faces = patch().localFaces();
    const labelListList& pointFaces = patch().pointFaces();

    for (const label curPoint : intPoints)
    {
        faceList curFaceList(pointFaces[curPoint].size());

        forAll(curFaceList, faceI)
        {
            curFaceList[faceI] = faces[pointFaces[curPoint][faceI]];
        }

        primitiveFacePatch curPatch(curFaceList, points);

        labelList curPointPoints = curPatch.edgeLoops()[0];

        for (label i = 0; i < curPointPoints.size(); ++i)
        {
            vector d1 =
                points[curPatch.meshPoints()[curPointPoints[i]]]
              - points[curPoint];

            label p = i + 1;

            if (i == (curPointPoints.size() - 1))
            {
                p = 0;
            }

            vector d2 =
                points[curPatch.meshPoints()[curPointPoints[p]]]
              - points[curPoint];

            vector n = (d1 ^ d2)/(mag(d1 ^ d2) + SMALL);

            scalar sinAlpha = mag(d1^d2)/(mag(d1)*mag(d2));

            scalar w = sinAlpha/(mag(d1)*mag(d2));

            result[curPoint] += w*n;
        }
    }

    for (const label curPoint : bndPoints)
    {
        faceList curFaceList(pointFaces[curPoint].size());

        forAll(curFaceList, faceI)
        {
            curFaceList[faceI] = faces[pointFaces[curPoint][faceI]];
        }

        primitiveFacePatch curPatch(curFaceList, points);

        labelList agglomFacePoints = curPatch.edgeLoops()[0];

        SLList<label> slList;

        label curPointLabel = -1;

        for (label i=0; i<agglomFacePoints.size(); ++i)
        {
            if (curPatch.meshPoints()[agglomFacePoints[i]] == curPoint)
            {
                curPointLabel = i;
            }
            else if ( curPointLabel != -1 )
            {
                slList.append(curPatch.meshPoints()[agglomFacePoints[i]]);
            }
        }

        for (label i=0; i<curPointLabel; ++i)
        {
            slList.append(curPatch.meshPoints()[agglomFacePoints[i]]);
        }

        labelList curPointPoints(slList);

        for (label i=0; i < (curPointPoints.size() - 1); ++i)
        {
            vector d1 = points[curPointPoints[i]] - points[curPoint];

            vector d2 = points[curPointPoints[i + 1]] - points[curPoint];

            vector n = (d1 ^ d2)/(mag(d1 ^ d2) + SMALL);

            scalar sinAlpha = mag(d1 ^ d2)/(mag(d1)*mag(d2));

            scalar w = sinAlpha/(mag(d1)*mag(d2));

            result[curPoint] += w*n;
        }
    }

    // Correct wedge points
    forAll(boundary(), patchI)
    {
        const faPatch& fap = boundary()[patchI];

        if (isA<wedgeFaPatch>(fap))
        {
            const wedgeFaPatch& wedgePatch = refCast<const wedgeFaPatch>(fap);

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
                result[pti] -= N*(N&result[pti]);
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
        }
    }


    // Processor patch points correction
    for (const faPatch& fap : boundary())
    {
        if (Pstream::parRun() && isA<processorFaPatch>(fap))
        {
            const processorFaPatch& procPatch =
                refCast<const processorFaPatch>(fap);

            const labelList& patchPointLabels = procPatch.pointLabels();

            vectorField patchPointNormals
            (
                patchPointLabels.size(),
                Zero
            );

            forAll(patchPointNormals, pointI)
            {
                patchPointNormals[pointI] = result[patchPointLabels[pointI]];
            }

            {
            OPstream::write
            (
                Pstream::commsTypes::blocking,
                procPatch.neighbProcNo(),
                patchPointNormals.cdata_bytes(),
                patchPointNormals.byteSize()
            );
            }

            vectorField ngbPatchPointNormals
            (
                procPatch.neighbPoints().size(),
                Zero
            );

            {
                IPstream::read
                (
                    Pstream::commsTypes::blocking,
                    procPatch.neighbProcNo(),
                    ngbPatchPointNormals.data_bytes(),
                    ngbPatchPointNormals.byteSize()
                );
            }

            const labelList& nonGlobalPatchPoints =
                procPatch.nonGlobalPatchPoints();

            for (const label pti : nonGlobalPatchPoints)
            {
                result[patchPointLabels[pti]] +=
                    ngbPatchPointNormals[procPatch.neighbPoints()[pti]];
            }
        }
    }


    // Correct global points
    if (globalData().nGlobalPoints() > 0)
    {
        const labelList& spLabels(globalData().sharedPointLabels());

        vectorField spNormals(spLabels.size(), Zero);
        forAll(spNormals, pointI)
        {
            spNormals[pointI] = result[spLabels[pointI]];
        }

        const labelList& addr = globalData().sharedPointAddr();

        vectorField gpNormals
        (
            globalData().nGlobalPoints(),
            Zero
        );

        forAll(addr, i)
        {
            gpNormals[addr[i]] += spNormals[i];
        }

        combineReduce(gpNormals, plusEqOp<vectorField>());

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


    // Boundary points correction
    forAll(boundary(), patchI)
    {
        const faPatch& fap = boundary()[patchI];

        if (correctPatchPointNormals(patchI) && !fap.coupled())
        {
            if (fap.ngbPolyPatchIndex() < 0)
            {
                FatalErrorInFunction
                    << "Neighbour polyPatch index is not defined "
                    << "for faPatch " << fap.name()
                    << abort(FatalError);
            }

            const labelList& patchPoints = fap.pointLabels();

            const vectorField N(fap.ngbPolyPatchPointNormals());

            forAll(patchPoints, pointI)
            {
                result[patchPoints[pointI]]
                    -= N[pointI]*(N[pointI]&result[patchPoints[pointI]]);
            }
        }
    }

    for (vector& n : result)
    {
        n.normalise();
    }
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

    result.resize(nPoints());
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
                result[pti] -= N*(N&result[pti]);
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
                    reinterpret_cast<const char*>(patchPointNormals.cdata()),
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
                    reinterpret_cast<char*>(patchPointNormals.data()),
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

        combineReduce(gpNormals, plusEqOp<vectorField>());

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
            vector fpnorm = iter.val();

            fpnorm.normalise();
            result[pointi] -= fpnorm*(fpnorm & result[pointi]);
        }
    }

    for (vector& n : result)
    {
        n.normalise();
    }
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
                    toNgbProcLsPoints.byteSize()
                  + toNgbProcLsPointStarts.byteSize()
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
                  + fromNgbProcLsPointStarts.byteSize()
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
                dir -= axis*(axis&dir);
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

            Pstream::gatherList(procLsPoints);
            Pstream::scatterList(procLsPoints);

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
                dir -= axis*(axis&dir);
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

    tmp<edgeScalarField> tcorrection
    (
        new edgeScalarField
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
        )
    );
    edgeScalarField& correction = tcorrection.ref();

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
