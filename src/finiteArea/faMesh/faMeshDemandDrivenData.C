/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2016-2017 Wikki Ltd
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
#include "coordinateSystem.H"
#include "scalarMatrices.H"
#include "processorFaPatchFields.H"
#include "emptyFaPatchFields.H"

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

    patchStartsPtr_ = new labelList(boundary().size(), -1);
    labelList& patchStarts = *patchStartsPtr_;

    patchStarts[0] = nInternalEdges();

    for (label i = 1; i < boundary().size(); ++i)
    {
        patchStarts[i] =
            patchStarts[i - 1] + boundary()[i - 1].faPatch::size();
    }
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

    const edgeList::subList internalEdges =
        edgeList::subList(edges(), nInternalEdges());


    forAll(internalEdges, edgeI)
    {
        magLe.ref()[edgeI] = internalEdges[edgeI].mag(localPoints);
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

    const edgeList::subList internalEdges =
        edgeList::subList(edges(), nInternalEdges());


    forAll(internalEdges, edgeI)
    {
        edgeCentres.ref()[edgeI] = internalEdges[edgeI].centre(localPoints);
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
        nInternal[faceI] =
            localFaces[faceI].normal(localPoints)/
            localFaces[faceI].mag(localPoints);
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


//     // Primitive patch edge normals
//     const labelListList& patchPointEdges = patch().pointEdges();

//     vectorField patchEdgeNormals(nEdges(), vector::zero);

//     forAll(pointNormals, pointI)
//     {
//         const labelList& curPointEdges = patchPointEdges[pointI];

//         forAll(curPointEdges, edgeI)
//         {
//             label curEdge = curPointEdges[edgeI];

//             patchEdgeNormals[curEdge] += 0.5*pointNormals[pointI];
//         }
//     }

//     patchEdgeNormals /= mag(patchEdgeNormals);


//     // Edge area normals
//     label nIntEdges = patch().nInternalEdges();

//     for (label edgeI = 0; edgeI < nIntEdges; ++edgeI)
//     {
//         edgeAreaNormals.ref()[edgeI] =
//             patchEdgeNormals[edgeI];
//     }

//     forAll(boundary(), patchI)
//     {
//         const labelList& edgeLabels = boundary()[patchI];

//         forAll(edgeAreaNormals.boundaryFieldRef()[patchI], edgeI)
//         {
//             edgeAreaNormals.boundaryFieldRef()[patchI][edgeI] =
//                 patchEdgeNormals[edgeLabels[edgeI]];
//         }
//     }


    forAll(edgeAreaNormals.internalField(), edgeI)
    {
        vector e = edges()[edgeI].vec(points());
        e /= mag(e);

//         scalar wStart =
//             1.0 - sqr(mag(e^pointNormals[edges()[edgeI].end()]));

//         scalar wEnd =
//             1.0 - sqr(mag(e^pointNormals[edges()[edgeI].start()]));

//         wStart = 1.0;
//         wEnd = 1.0;

//         edgeAreaNormals.ref()[edgeI] =
//             wStart*pointNormals[edges()[edgeI].start()]
//           + wEnd*pointNormals[edges()[edgeI].end()];

//         vector eC =
//             0.5
//            *(
//                  points()[edges()[edgeI].start()]
//                + points()[edges()[edgeI].end()]
//             );

//         vector eCp = 0.5*
//             (
//                 points()[edges()[edgeI].start()]
//               + pointNormals[edges()[edgeI].start()]
//                 points()[edges()[edgeI].end()] +
//             );

        edgeAreaNormals.ref()[edgeI] =
            pointNormals[edges()[edgeI].start()]
          + pointNormals[edges()[edgeI].end()];

        edgeAreaNormals.ref()[edgeI] -=
            e*(e&edgeAreaNormals.internalField()[edgeI]);
    }

    edgeAreaNormals.ref() /= mag(edgeAreaNormals.internalField());

    forAll(boundary(), patchI)
    {
        const edgeList::subList patchEdges =
            boundary()[patchI].patchSlice(edges());

        forAll(patchEdges, edgeI)
        {
            edgeAreaNormals.boundaryFieldRef()[patchI][edgeI] =
                pointNormals[patchEdges[edgeI].start()]
              + pointNormals[patchEdges[edgeI].end()];

            vector e = patchEdges[edgeI].vec(points());
            e /= mag(e);

            edgeAreaNormals.boundaryFieldRef()[patchI][edgeI] -=
                e*(e&edgeAreaNormals.boundaryField()[patchI][edgeI]);
        }

        edgeAreaNormals.boundaryFieldRef()[patchI] /=
            mag(edgeAreaNormals.boundaryField()[patchI]);
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

    const edgeList& edges = patch().edges();
    label nIntEdges = patch().nInternalEdges();

    List<bool> internal(nPoints(), true);

    for (label curEdge = nIntEdges; curEdge < edges.size(); ++curEdge)
    {
        internal[edges[curEdge].start()] = false;

        internal[edges[curEdge].end()] = false;
    }

    SLList<label> internalPoints;

    forAll(internal, pointI)
    {
        if (internal[pointI])
        {
            internalPoints.append(pointI);
        }
    }

    labelList result(internalPoints);

    return result;
}


Foam::labelList Foam::faMesh::boundaryPoints() const
{
    DebugInFunction
        << "Calculating boundary points" << endl;

    const edgeList& edges = patch().edges();
    label nIntEdges = patch().nInternalEdges();

    List<bool> internal(nPoints(), true);

    for (label curEdge = nIntEdges; curEdge < edges.size(); ++curEdge)
    {
        internal[edges[curEdge].start()] = false;

        internal[edges[curEdge].end()] = false;
    }

    SLList<label> boundaryPoints;

    forAll(internal, pointI)
    {
        if (!internal[pointI])
        {
            boundaryPoints.append(pointI);
        }
    }

    labelList result(boundaryPoints);

    return result;
}


void Foam::faMesh::calcPointAreaNormals() const
{
    if (pointAreaNormalsPtr_)
    {
        FatalErrorInFunction
            << "pointAreaNormalsPtr_ already allocated"
            << abort(FatalError);
    }


    pointAreaNormalsPtr_ = new vectorField(nPoints(), vector::zero);

    vectorField& result = *pointAreaNormalsPtr_;

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

        for (int i = 0; i < curPointPoints.size(); ++i)
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

            vector N =
                transform
                (
                    wedgePatch.edgeT(),
                    wedgePatch.centreNormal()
                );

            N /= mag(N);

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
        if (isA<processorFaPatch>(fap))
        {
            const processorFaPatch& procPatch =
                refCast<const processorFaPatch>(fap);

            const labelList& patchPointLabels = procPatch.pointLabels();

            vectorField patchPointNormals
            (
                patchPointLabels.size(),
                vector::zero
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
                reinterpret_cast<const char*>(patchPointNormals.begin()),
                patchPointNormals.byteSize()
            );
            }

            vectorField ngbPatchPointNormals
            (
                procPatch.neighbPoints().size(),
                vector::zero
            );

            {
                IPstream::read
                (
                    Pstream::commsTypes::blocking,
                    procPatch.neighbProcNo(),
                    reinterpret_cast<char*>(ngbPatchPointNormals.begin()),
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

        vectorField spNormals(spLabels.size(), vector::zero);
        forAll(spNormals, pointI)
        {
            spNormals[pointI] = result[spLabels[pointI]];
        }

        const labelList& addr = globalData().sharedPointAddr();

        vectorField gpNormals
        (
            globalData().nGlobalPoints(),
            vector::zero
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
            if (fap.ngbPolyPatchIndex() == -1)
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

    result /= mag(result);
}


void Foam::faMesh::calcPointAreaNormalsByQuadricsFit() const
{
    vectorField& result = *pointAreaNormalsPtr_;

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

        pointSet.insert(curPoint);
        for (label i=0; i<curFaces.size(); ++i)
        {
            const labelList& facePoints = faces[curFaces[i]];
            for (label j=0; j<facePoints.size(); ++j)
            {
                if (!pointSet.found(facePoints[j]))
                {
                    pointSet.insert(facePoints[j]);
                }
            }
        }
        pointSet.erase(curPoint);
        labelList curPoints(pointSet.toc());

        if (curPoints.size() < 5)
        {
            if (debug)
            {
                Info << "WARNING: Extending point set for fitting." << endl;
            }

            labelHashSet faceSet(pointFaces[curPoint]);
            labelList curFaces(faceSet.toc());
            forAll(curFaces, faceI)
            {
                const labelList& curFaceFaces =
                    patch().faceFaces()[curFaces[faceI]];

                forAll(curFaceFaces, fI)
                {
                    label curFaceFace = curFaceFaces[fI];

                    if (!faceSet.found(curFaceFace))
                    {
                        faceSet.insert(curFaceFace);
                    }
                }
            }
            curFaces = faceSet.toc();

            labelHashSet pointSet;

            pointSet.insert(curPoint);
            for (label i=0; i<curFaces.size(); ++i)
            {
                const labelList& facePoints = faces[curFaces[i]];
                for (label j=0; j<facePoints.size(); ++j)
                {
                    if (!pointSet.found(facePoints[j]))
                    {
                        pointSet.insert(facePoints[j]);
                    }
                }
            }

            pointSet.erase(curPoint);
            curPoints = pointSet.toc();
        }

        vectorField allPoints(curPoints.size());
        scalarField W(curPoints.size(), 1.0);
        for (label i=0; i<curPoints.size(); ++i)
        {
            allPoints[i] = points[curPoints[i]];
            W[i] = 1.0/magSqr(allPoints[i] - points[curPoint]);
        }

        // Transforme points
        const vector& origin = points[curPoint];
        vector axis(result[curPoint]/mag(result[curPoint]));
        vector dir(allPoints[0] - points[curPoint]);
        dir -= axis*(axis&dir);
        dir /= mag(dir);
        coordinateSystem cs("cs", origin, axis, dir);

        forAll(allPoints, pI)
        {
            allPoints[pI] = cs.localPosition(allPoints[pI]);
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

        scalarSquareMatrix MtM(5, 0.0);

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

        scalarField MtR(5, 0);

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

        if (isA<processorFaPatch>(fap))
        {
            const processorFaPatch& procPatch =
                refCast<const processorFaPatch>(boundary()[patchI]);

            const labelList& patchPointLabels = procPatch.pointLabels();

            labelList toNgbProcLsPointStarts(patchPointLabels.size(), 0);
            vectorField toNgbProcLsPoints
            (
                10*patchPointLabels.size(),
                vector::zero
            );
            label nPoints = 0;

            for (label pointI=0; pointI<patchPointLabels.size(); ++pointI)
            {
                label curPoint = patchPointLabels[pointI];

                toNgbProcLsPointStarts[pointI] = nPoints;

                const labelHashSet faceSet(pointFaces[curPoint]);
                const labelList curFaces(faceSet.toc());

                labelHashSet pointSet;

                pointSet.insert(curPoint);
                for (label i=0; i<curFaces.size(); ++i)
                {
                    const labelList& facePoints = faces[curFaces[i]];
                    for (label j=0; j<facePoints.size(); ++j)
                    {
                        if (!pointSet.found(facePoints[j]))
                        {
                            pointSet.insert(facePoints[j]);
                        }
                    }
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
        if (isA<processorFaPatch>(fap))
        {
            const processorFaPatch& procPatch =
                refCast<const processorFaPatch>(fap);

            const labelList& patchPointLabels = procPatch.pointLabels();

            labelList fromNgbProcLsPointStarts(patchPointLabels.size(), 0);
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
                forAll(pointFaces[curPoint], faceI)
                {
                    faceSet.insert(pointFaces[curPoint][faceI]);
                }
                labelList curFaces = faceSet.toc();

                labelHashSet pointSet;

                pointSet.insert(curPoint);
                for (label i=0; i<curFaces.size(); ++i)
                {
                    const labelList& facePoints = faces[curFaces[i]];
                    for (label j=0; j<facePoints.size(); ++j)
                    {
                        if (!pointSet.found(facePoints[j]))
                        {
                            pointSet.insert(facePoints[j]);
                        }
                    }
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
                vectorField allPoints(nAllPoints, vector::zero);
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
                vector axis(result[curPoint]/mag(result[curPoint]));
                vector dir(allPoints[0] - points[curPoint]);
                dir -= axis*(axis&dir);
                dir /= mag(dir);
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

                scalarSquareMatrix MtM(5, 0.0);

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

                scalarField MtR(5, 0);

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

            label curSharedPointIndex = findIndex(addr, k);

            scalar tol = 0.0;

            if (curSharedPointIndex != -1)
            {
                label curPoint = spLabels[curSharedPointIndex];

                const labelHashSet faceSet(pointFaces[curPoint]);
                const labelList curFaces(faceSet.toc());

                labelHashSet pointSet;
                pointSet.insert(curPoint);
                for (label i=0; i<curFaces.size(); ++i)
                {
                    const labelList& facePoints = faces[curFaces[i]];
                    for (label j=0; j<facePoints.size(); ++j)
                    {
                        if (!pointSet.found(facePoints[j]))
                        {
                            pointSet.insert(facePoints[j]);
                        }
                    }
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

                vectorField allPoints(nAllPoints, vector::zero);

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
                vector axis(result[curPoint]/mag(result[curPoint]));
                vector dir(allPoints[0] - points[curPoint]);
                dir -= axis*(axis&dir);
                dir /= mag(dir);
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

                scalarSquareMatrix MtM(5, 0.0);
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

                scalarField MtR(5, 0);
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

    result /= mag(result);
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
