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

#include "cyclicPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "demandDrivenData.H"
#include "OFstream.H"
#include "patchZones.H"
#include "matchPoints.H"
#include "EdgeMap.H"
#include "Time.H"
#include "diagTensor.H"
#include "transformField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, cyclicPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, cyclicPolyPatch, dictionary);


template<>
const char* NamedEnum<cyclicPolyPatch::transformType, 3>::names[] =
{
    "unknown",
    "rotational",
    "translational"
};

const NamedEnum<cyclicPolyPatch::transformType, 3>
    cyclicPolyPatch::transformTypeNames;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::label Foam::cyclicPolyPatch::findMaxArea
(
    const pointField& points,
    const faceList& faces
)
{
    label maxI = -1;
    scalar maxAreaSqr = -GREAT;

    forAll(faces, faceI)
    {
        scalar areaSqr = magSqr(faces[faceI].normal(points));

        if (areaSqr > maxAreaSqr)
        {
            maxAreaSqr = areaSqr;
            maxI = faceI;
        }
    }
    return maxI;
}


void Foam::cyclicPolyPatch::calcTransforms()
{
    if (size())
    {
        // Half0

        const cyclicPolyPatch& half0 = *this;

        const pointField& half0Ctrs = half0.faceCentres();

        if (debug)
        {
            fileName casePath(boundaryMesh().mesh().time().path());

            fileName nm0(casePath/name()+"_half0_faces.obj");
            Pout<< "cyclicPolyPatch::calcTransforms : Writing half0"
                << " faces to OBJ file " << nm0 << endl;
            writeOBJ(nm0, half0, half0.points());
        }

        vectorField half0Areas(half0.size());

        forAll(half0, facei)
        {
            half0Areas[facei] = half0[facei].normal(half0.points());
        }



        // Half0

        const cyclicPolyPatch& half1 = neighbPatch();

        const pointField& half1Ctrs = half1.faceCentres();

        // Dump halves
        if (debug)
        {
            fileName casePath(boundaryMesh().mesh().time().path());

            fileName nm1(casePath/name()+"_half1_faces.obj");
            Pout<< "cyclicPolyPatch::calcTransforms : Writing half1"
                << " faces to OBJ file " << nm1 << endl;
            writeOBJ(nm1, half1, half1.points());

            OFstream str(casePath/name()+"_half0_to_half1.obj");
            label vertI = 0;
            Pout<< "cyclicPolyPatch::calcTransforms :"
                << " Writing coupled face centres as lines to " << str.name()
                << endl;
            forAll(half0Ctrs, i)
            {
                const point& p0 = half0Ctrs[i];
                str << "v " << p0.x() << ' ' << p0.y() << ' ' << p0.z() << nl;
                vertI++;
                const point& p1 = half1Ctrs[i];
                str << "v " << p1.x() << ' ' << p1.y() << ' ' << p1.z() << nl;
                vertI++;
                str << "l " << vertI-1 << ' ' << vertI << nl;
            }
        }

        vectorField half1Areas(half1.size());

        forAll(half1, facei)
        {
            half1Areas[facei] = half1[facei].normal(half1.points());
        }

        calcTransforms
        (
            half0,
            half0Ctrs,
            half0Areas,
            half1Ctrs,
            half1Areas
        );
    }
}


void Foam::cyclicPolyPatch::calcTransforms
(
    const primitivePatch& half0,
    const UList<point>& half0Ctrs,
    const UList<point>& half0Areas,
    const UList<point>& half1Ctrs,
    const UList<point>& half1Areas
)
{
Pout<< "cyclicPolyPatch::calcTransforms : name:" << name() << endl
    << "    half0Ctrs:"
    << " min:" << min(half0Ctrs) << " max:" << max(half0Ctrs)<< endl
    << "    half1Ctrs:"
    << " min:" << min(half1Ctrs) << " max:" << max(half1Ctrs)<< endl
    << endl;

    if (half0Ctrs.size() > 0)
    {
        scalarField half0Tols
        (
            calcFaceTol
            (
                half0,
                half0.points(),
                static_cast<const pointField&>(half0Ctrs)
            )
        );

        vectorField half0Normals(half0Areas.size());
        vectorField half1Normals(half1Areas.size());

        forAll(half0, facei)
        {
            scalar magSf = mag(half0Areas[facei]);
            scalar nbrMagSf = mag(half1Areas[facei]);
            scalar avSf = (magSf + nbrMagSf)/2.0;

            if (magSf < ROOTVSMALL && nbrMagSf < ROOTVSMALL)
            {
                // Undetermined normal. Use dummy normal to force separation
                // check. (note use of sqrt(VSMALL) since that is how mag
                // scales)
                half0Normals[facei] = point(1, 0, 0);
                half1Normals[facei] = half0Normals[facei];
            }
            else if (mag(magSf - nbrMagSf)/avSf > coupledPolyPatch::matchTol)
            {
                FatalErrorIn
                (
                    "cyclicPolyPatch::calcTransforms()"
                )   << "face " << facei << " area does not match neighbour by "
                    << 100*mag(magSf - nbrMagSf)/avSf
                    << "% -- possible face ordering problem." << endl
                    << "patch:" << name()
                    << " my area:" << magSf
                    << " neighbour area:" << nbrMagSf
                    << " matching tolerance:" << coupledPolyPatch::matchTol
                     << endl
                    << "Mesh face:" << start()+facei
                    << " fc:" << half0Ctrs[facei]
                    << endl
                    << "Neighbour fc:" << half1Ctrs[facei]
                    << endl
                    << "Rerun with cyclic debug flag set"
                    << " for more information." << exit(FatalError);
            }
            else
            {
                half0Normals[facei] = half0Areas[facei] / magSf;
                half1Normals[facei] = half1Areas[facei] / nbrMagSf;
            }
        }

        // Calculate transformation tensors
        calcTransformTensors
        (
            separated_,
            separation_,
            parallel_,
            forwardT_,
            reverseT_,
            static_cast<const pointField&>(half0Ctrs),
            static_cast<const pointField&>(half1Ctrs),
            half0Normals,
            half1Normals,
            half0Tols
        );

Pout<< "cyclicPolyPatch::calcTransforms : calculated transforms for:"
    << name() << endl
    << "    separated_:" << separated_ << endl
    << "    separation_:" << separation_ << endl
    << "    parallel_:" << parallel_ << endl
    << "    forwardT_:" << forwardT_ << endl
    << "    reverseT_:" << reverseT_ << endl
    << endl;

    }
}


// Given a split of faces into left and right half calculate the centres
// and anchor points. Transform the left points so they align with the
// right ones.
void Foam::cyclicPolyPatch::getCentresAndAnchors
(
    const primitivePatch& pp0,
    const primitivePatch& pp1,

    pointField& half0Ctrs,
    pointField& half1Ctrs,
    pointField& anchors0,
    scalarField& tols
) const
{
    // Get geometric data on both halves.
    half0Ctrs = pp0.faceCentres();
    anchors0 = getAnchorPoints(pp0, pp0.points());
    half1Ctrs = pp1.faceCentres();

    switch (transform_)
    {
        case ROTATIONAL:
        {
            label face0 = getConsistentRotationFace(half0Ctrs);
            label face1 = getConsistentRotationFace(half1Ctrs);

            vector n0 = ((half0Ctrs[face0] - rotationCentre_) ^ rotationAxis_);
            vector n1 = ((half1Ctrs[face1] - rotationCentre_) ^ -rotationAxis_);
            n0 /= mag(n0) + VSMALL;
            n1 /= mag(n1) + VSMALL;

            if (debug)
            {
                Pout<< "cyclicPolyPatch::getCentresAndAnchors :"
                    << " Specified rotation :"
                    << " n0:" << n0 << " n1:" << n1 << endl;
            }

            // Rotation (around origin)
            const tensor reverseT(rotationTensor(n0, -n1));

            // Rotation
            forAll(half0Ctrs, faceI)
            {
                half0Ctrs[faceI] = Foam::transform(reverseT, half0Ctrs[faceI]);
                anchors0[faceI] = Foam::transform(reverseT, anchors0[faceI]);
            }

            break;
        }
        //- Problem: usually specified translation is not accurate enough
        //- to get proper match so keep automatic determination over here.
        //case TRANSLATIONAL:
        //{
        //    // Transform 0 points.
        //
        //    if (debug)
        //    {
        //        Pout<< "cyclicPolyPatch::getCentresAndAnchors :"
        //            << "Specified translation : " << separationVector_
        //            << endl;
        //    }
        //
        //    half0Ctrs += separationVector_;
        //    anchors0 += separationVector_;
        //    break;
        //}
        default:
        {
            // Assumes that cyclic is planar. This is also the initial
            // condition for patches without faces.

            // Determine the face with max area on both halves. These
            // two faces are used to determine the transformation tensors
            label max0I = findMaxArea(pp0.points(), pp0);
            vector n0 = pp0[max0I].normal(pp0.points());
            n0 /= mag(n0) + VSMALL;

            label max1I = findMaxArea(pp1.points(), pp1);
            vector n1 = pp1[max1I].normal(pp1.points());
            n1 /= mag(n1) + VSMALL;

            if (mag(n0 & n1) < 1-coupledPolyPatch::matchTol)
            {
                if (debug)
                {
                    Pout<< "cyclicPolyPatch::getCentresAndAnchors :"
                        << " Detected rotation :"
                        << " n0:" << n0 << " n1:" << n1 << endl;
                }

                // Rotation (around origin)
                const tensor reverseT(rotationTensor(n0, -n1));

                // Rotation
                forAll(half0Ctrs, faceI)
                {
                    half0Ctrs[faceI] = Foam::transform
                    (
                        reverseT,
                        half0Ctrs[faceI]
                    );
                    anchors0[faceI] = Foam::transform
                    (
                        reverseT,
                        anchors0[faceI]
                    );
                }
            }
            else
            {
                // Parallel translation. Get average of all used points.

                const point ctr0(sum(pp0.localPoints())/pp0.nPoints());
                const point ctr1(sum(pp1.localPoints())/pp1.nPoints());

                if (debug)
                {
                    Pout<< "cyclicPolyPatch::getCentresAndAnchors :"
                        << " Detected translation :"
                        << " n0:" << n0 << " n1:" << n1
                        << " ctr0:" << ctr0 << " ctr1:" << ctr1 << endl;
                }

                half0Ctrs += ctr1 - ctr0;
                anchors0 += ctr1 - ctr0;
            }
            break;
        }
    }


    // Calculate typical distance per face
    tols = calcFaceTol(pp1, pp1.points(), half1Ctrs);
}


Foam::label Foam::cyclicPolyPatch::getConsistentRotationFace
(
    const pointField& faceCentres
) const
{
    const scalarField magRadSqr =
        magSqr((faceCentres - rotationCentre_) ^ rotationAxis_);
    scalarField axisLen = (faceCentres - rotationCentre_) & rotationAxis_;
    axisLen = axisLen - min(axisLen);
    const scalarField magLenSqr = magRadSqr + axisLen*axisLen;

    label rotFace = -1;
    scalar maxMagLenSqr = -GREAT;
    scalar maxMagRadSqr = -GREAT;
    forAll(faceCentres, i)
    {
        if (magLenSqr[i] >= maxMagLenSqr)
        {
            if (magRadSqr[i] > maxMagRadSqr)
            {
                rotFace = i;
                maxMagLenSqr = magLenSqr[i];
                maxMagRadSqr = magRadSqr[i];
            }
        }
    }

    if (debug)
    {
        Info<< "getConsistentRotationFace(const pointField&)" << nl
            << "    rotFace = " << rotFace << nl
            << "    point =  " << faceCentres[rotFace] << endl;
    }

    return rotFace;
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::cyclicPolyPatch::cyclicPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, size, start, index, bm),
    neighbPatchName_(word::null),
    neighbPatchID_(-1),
    transform_(UNKNOWN),
    rotationAxis_(vector::zero),
    rotationCentre_(point::zero),
    separationVector_(vector::zero),
    coupledPointsPtr_(NULL),
    coupledEdgesPtr_(NULL),
    separated_(false),
    separation_(vector::zero),
    parallel_(true),
    forwardT_(I),
    reverseT_(I)
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible.
}


Foam::cyclicPolyPatch::cyclicPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, dict, index, bm),
    neighbPatchName_(dict.lookup("neighbourPatch")),
    neighbPatchID_(-1),
    transform_(UNKNOWN),
    rotationAxis_(vector::zero),
    rotationCentre_(point::zero),
    separationVector_(vector::zero),
    coupledPointsPtr_(NULL),
    coupledEdgesPtr_(NULL),
    separated_(false),
    separation_(vector::zero),
    parallel_(true),
    forwardT_(I),
    reverseT_(I)
{
    if (dict.found("transform"))
    {
        transform_ = transformTypeNames.read(dict.lookup("transform"));
        switch (transform_)
        {
            case ROTATIONAL:
            {
                dict.lookup("rotationAxis") >> rotationAxis_;
                dict.lookup("rotationCentre") >> rotationCentre_;
                break;
            }
            case TRANSLATIONAL:
            {
                dict.lookup("separationVector") >> separationVector_;
                break;
            }
            default:
            {
                // no additional info required
            }
        }
    }

    // Neighbour patch might not be valid yet so no transformation
    // calculation possible.
}


Foam::cyclicPolyPatch::cyclicPolyPatch
(
    const cyclicPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    neighbPatchName_(pp.neighbPatchName()),
    neighbPatchID_(-1),
    transform_(pp.transform_),
    rotationAxis_(pp.rotationAxis_),
    rotationCentre_(pp.rotationCentre_),
    separationVector_(pp.separationVector_),
    coupledPointsPtr_(NULL),
    coupledEdgesPtr_(NULL),
    separated_(false),
    separation_(vector::zero),
    parallel_(true),
    forwardT_(I),
    reverseT_(I)
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible.
}


Foam::cyclicPolyPatch::cyclicPolyPatch
(
    const cyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& neighbPatchName
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    neighbPatchName_(neighbPatchName),
    neighbPatchID_(-1),
    transform_(pp.transform_),
    rotationAxis_(pp.rotationAxis_),
    rotationCentre_(pp.rotationCentre_),
    separationVector_(pp.separationVector_),
    coupledPointsPtr_(NULL),
    coupledEdgesPtr_(NULL)
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible.
}


Foam::cyclicPolyPatch::cyclicPolyPatch
(
    const cyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const unallocLabelList& mapAddressing,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, mapAddressing, newStart),
    neighbPatchName_(pp.neighbPatchName_),
    neighbPatchID_(-1),
    transform_(pp.transform_),
    rotationAxis_(pp.rotationAxis_),
    rotationCentre_(pp.rotationCentre_),
    separationVector_(pp.separationVector_),
    coupledPointsPtr_(NULL),
    coupledEdgesPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicPolyPatch::~cyclicPolyPatch()
{
    deleteDemandDrivenData(coupledPointsPtr_);
    deleteDemandDrivenData(coupledEdgesPtr_);
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cyclicPolyPatch::transformPosition(pointField& l) const
{
    if (!parallel())
    {
        Foam::transform(forwardT_, l);
    }
    else if (separated())
    {
        l -= separation_;
    }
}


void Foam::cyclicPolyPatch::initGeometry(PstreamBuffers& pBufs)
{
    polyPatch::initGeometry(pBufs);
}


void Foam::cyclicPolyPatch::initGeometry
(
    const primitivePatch& referPatch,
    UList<point>& nbrCtrs,
    UList<point>& nbrAreas,
    UList<point>& nbrCc
)
{}


void Foam::cyclicPolyPatch::calcGeometry
(
    const primitivePatch& referPatch,
    const UList<point>& thisCtrs,
    const UList<point>& thisAreas,
    const UList<point>& thisCc,
    const UList<point>& nbrCtrs,
    const UList<point>& nbrAreas,
    const UList<point>& nbrCc
)
{
    //polyPatch::calcGeometry();

Pout<< "cyclicPolyPatch::calcGeometry : name:" << name()
    << " referred from:" << referPatch.size() << endl;

    calcTransforms
    (
        referPatch,
        thisCtrs,
        thisAreas,
        nbrCtrs,
        nbrAreas
    );
}


void Foam::cyclicPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::initMovePoints(pBufs, p);
}

void Foam::cyclicPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::movePoints(pBufs, p);
    calcTransforms();
}

void Foam::cyclicPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    polyPatch::initUpdateMesh(pBufs);
}

void Foam::cyclicPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    polyPatch::updateMesh(pBufs);
    deleteDemandDrivenData(coupledPointsPtr_);
    deleteDemandDrivenData(coupledEdgesPtr_);
}


const Foam::edgeList& Foam::cyclicPolyPatch::coupledPoints() const
{
    if (!coupledPointsPtr_)
    {
        const faceList& nbrLocalFaces = neighbPatch().localFaces();
        const labelList& nbrMeshPoints = neighbPatch().meshPoints();

        // Now all we know is that relative face index in *this is same
        // as coupled face in nbrPatch and also that the 0th vertex
        // corresponds.

        // From local point to nbrPatch or -1.
        labelList coupledPoint(nPoints(), -1);

        forAll(*this, patchFaceI)
        {
            const face& fA = localFaces()[patchFaceI];

            forAll(fA, indexA)
            {
                label patchPointA = fA[indexA];

                if (coupledPoint[patchPointA] == -1)
                {
                    const face& fB = nbrLocalFaces[patchFaceI];

                    label indexB = (fB.size() - indexA) % fB.size();

                    // Filter out points on wedge axis
                    if (patchPointA != fB[indexB])
                    {
                        coupledPoint[patchPointA] = fB[indexB];
                    }
                }
            }
        }

        coupledPointsPtr_ = new edgeList(nPoints());
        edgeList& connected = *coupledPointsPtr_;

        // Extract coupled points.
        label connectedI = 0;

        forAll(coupledPoint, i)
        {
            if (coupledPoint[i] != -1)
            {
                connected[connectedI++] = edge(i, coupledPoint[i]);
            }
        }

        connected.setSize(connectedI);

        if (debug)
        {
            OFstream str
            (
                boundaryMesh().mesh().time().path()
               /name() + "_coupledPoints.obj"
            );
            label vertI = 0;

            Pout<< "Writing file " << str.name() << " with coordinates of "
                << "coupled points" << endl;

            forAll(connected, i)
            {
                const point& a = points()[meshPoints()[connected[i][0]]];
                const point& b = points()[nbrMeshPoints[connected[i][1]]];

                str<< "v " << a.x() << ' ' << a.y() << ' ' << a.z() << nl;
                str<< "v " << b.x() << ' ' << b.y() << ' ' << b.z() << nl;
                vertI += 2;

                str<< "l " << vertI-1 << ' ' << vertI << nl;
            }
        }

        // Remove any addressing calculated for the coupled edges calculation
        const_cast<primitivePatch&>
        (
            static_cast<const primitivePatch&>
            (
                *this
            )
        ).clearOut();
        const_cast<primitivePatch&>
        (
            static_cast<const primitivePatch&>
            (
                neighbPatch()
            )
        ).clearOut();
    }
    return *coupledPointsPtr_;
}


const Foam::edgeList& Foam::cyclicPolyPatch::coupledEdges() const
{
    if (!coupledEdgesPtr_)
    {
        const edgeList& pointCouples = coupledPoints();

        // Build map from points on *this (A) to points on neighbourpatch (B)
        Map<label> aToB(2*pointCouples.size());

        forAll(pointCouples, i)
        {
            const edge& e = pointCouples[i];

            aToB.insert(e[0], e[1]);
        }

        // Map from edge on A to points (in B indices)
        EdgeMap<label> edgeMap(nEdges());

        forAll(*this, patchFaceI)
        {
            const labelList& fEdges = faceEdges()[patchFaceI];

            forAll(fEdges, i)
            {
                label edgeI = fEdges[i];

                const edge& e = edges()[edgeI];

                // Convert edge end points to corresponding points on B side.
                Map<label>::const_iterator fnd0 = aToB.find(e[0]);
                if (fnd0 != aToB.end())
                {
                    Map<label>::const_iterator fnd1 = aToB.find(e[1]);
                    if (fnd1 != aToB.end())
                    {
                        edgeMap.insert(edge(fnd0(), fnd1()), edgeI);
                    }
                }
            }
        }

        // Use the edgeMap to get the edges on the B side.

        const cyclicPolyPatch& neighbPatch = this->neighbPatch();

        coupledEdgesPtr_ = new edgeList(edgeMap.size());
        edgeList& coupledEdges = *coupledEdgesPtr_;
        label coupleI = 0;

        forAll(neighbPatch, patchFaceI)
        {
            const labelList& fEdges = neighbPatch.faceEdges()[patchFaceI];

            forAll(fEdges, i)
            {
                label edgeI = fEdges[i];

                const edge& e = neighbPatch.edges()[edgeI];

                // Look up A edge from HashTable.
                EdgeMap<label>::iterator iter = edgeMap.find(e);

                if (iter != edgeMap.end())
                {
                    label edgeA = iter();

                    // Store correspondence. Filter out edges on wedge axis.
                    if (edgeA != edgeI)
                    {
                        coupledEdges[coupleI++] = edge(edgeA, edgeI);
                    }

                    // Remove so we build unique list
                    edgeMap.erase(iter);
                }
            }
        }
        coupledEdges.setSize(coupleI);


        // Some checks

        forAll(coupledEdges, i)
        {
            const edge& e = coupledEdges[i];

            if (e[0] == e[1] || e[0] < 0 || e[1] < 0)
            {
                FatalErrorIn("cyclicPolyPatch::coupledEdges() const")
                    << "Problem : at position " << i
                    << " illegal couple:" << e
                    << abort(FatalError);
            }
        }

        if (debug)
        {
            OFstream str
            (
                boundaryMesh().mesh().time().path()
               /name() + "_coupledEdges.obj"
            );
            label vertI = 0;

            Pout<< "Writing file " << str.name() << " with centres of "
                << "coupled edges" << endl;

            forAll(coupledEdges, i)
            {
                const edge& e = coupledEdges[i];

                const point& a = edges()[e[0]].centre(localPoints());
                const point& b = neighbPatch.edges()[e[1]].centre
                (
                    neighbPatch.localPoints()
                );

                str<< "v " << a.x() << ' ' << a.y() << ' ' << a.z() << nl;
                str<< "v " << b.x() << ' ' << b.y() << ' ' << b.z() << nl;
                vertI += 2;

                str<< "l " << vertI-1 << ' ' << vertI << nl;
            }
        }

        // Remove any addressing calculated for the coupled edges calculation
        const_cast<primitivePatch&>
        (
            static_cast<const primitivePatch&>
            (
                *this
            )
        ).clearOut();
        const_cast<primitivePatch&>
        (
            static_cast<const primitivePatch&>
            (
                neighbPatch
            )
        ).clearOut();
    }
    return *coupledEdgesPtr_;
}


void Foam::cyclicPolyPatch::initOrder
(
    PstreamBuffers&,
    const primitivePatch& pp
) const
{
    if (owner())
    {
        // Save patch for use in non-owner side ordering. Equivalent to
        // processorPolyPatch using OPstream.
        ownerPatchPtr_.reset
        (
            new primitivePatch
            (
                pp,
                pp.points()
            )
        );
    }
}


//  Return new ordering. Ordering is -faceMap: for every face index
//  the new face -rotation:for every new face the clockwise shift
//  of the original face. Return false if nothing changes (faceMap
//  is identity, rotation is 0)
bool Foam::cyclicPolyPatch::order
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    faceMap.setSize(pp.size());
    faceMap = -1;

    rotation.setSize(pp.size());
    rotation = 0;

    if (pp.empty())
    {
        // No faces, nothing to change.
        return false;
    }

    if (owner())
    {
        // Do nothing (i.e. identical mapping, zero rotation).
        // See comment at top.
        forAll(faceMap, patchFaceI)
        {
            faceMap[patchFaceI] = patchFaceI;
        }

        return false;
    }
    else
    {
        // Get stored geometry from initOrder invocation of owner.
        const primitivePatch& pp0 = ownerPatchPtr_();

        // Get geometric quantities
        pointField half0Ctrs, half1Ctrs, anchors0;
        scalarField tols;
        getCentresAndAnchors
        (
            pp0,
            pp,

            half0Ctrs,
            half1Ctrs,
            anchors0,
            tols
        );

        // Geometric match of face centre vectors
        bool matchedAll = matchPoints
        (
            half1Ctrs,
            half0Ctrs,
            tols,
            true,
            faceMap
        );

        if (!matchedAll || debug)
        {
            // Dump halves
            fileName nm0
            (
                boundaryMesh().mesh().time().path()
               /name()+"_half0_faces.obj"
            );
            Pout<< "cyclicPolyPatch::order : Writing half0"
                << " faces to OBJ file " << nm0 << endl;
            writeOBJ(nm0, pp0, pp0.points());

            fileName nm1
            (
                boundaryMesh().mesh().time().path()
               /name()+"_half1_faces.obj"
            );
            Pout<< "cyclicPolyPatch::order : Writing half1"
                << " faces to OBJ file " << nm1 << endl;
            writeOBJ(nm1, pp, pp.points());

            OFstream ccStr
            (
                boundaryMesh().mesh().time().path()
               /name() + "_faceCentres.obj"
            );
            Pout<< "cyclicPolyPatch::order : "
                << "Dumping currently found cyclic match as lines between"
                << " corresponding face centres to file " << ccStr.name()
                << endl;

            // Recalculate untransformed face centres
            //pointField rawHalf0Ctrs =
            //    calcFaceCentres(half0Faces, pp.points());
            label vertI = 0;

            forAll(half1Ctrs, i)
            {
                if (faceMap[i] != -1)
                {
                    // Write edge between c1 and c0
                    const point& c0 = half0Ctrs[faceMap[i]];
                    const point& c1 = half1Ctrs[i];
                    writeOBJ(ccStr, c0, c1, vertI);
                }
            }
        }

        if (!matchedAll)
        {
            SeriousErrorIn
            (
                "cyclicPolyPatch::order"
                "(const primitivePatch&, labelList&, labelList&) const"
            )   << "Patch:" << name() << " : "
                << "Cannot match vectors to faces on both sides of patch"
                << endl
                << "    Perhaps your faces do not match?"
                << " The obj files written contain the current match." << endl
                << "    Continuing with incorrect face ordering from now on!"
                << endl;

                return false;
        }


        // Set rotation.
        forAll(faceMap, oldFaceI)
        {
            // The face f will be at newFaceI (after morphing) and we want its
            // anchorPoint (= f[0]) to align with the anchorpoint for the
            // corresponding face on the other side.

            label newFaceI = faceMap[oldFaceI];

            const point& wantedAnchor = anchors0[newFaceI];

            rotation[newFaceI] = getRotation
            (
                pp.points(),
                pp[oldFaceI],
                wantedAnchor,
                tols[oldFaceI]
            );

            if (rotation[newFaceI] == -1)
            {
                SeriousErrorIn
                (
                    "cyclicPolyPatch::order(const primitivePatch&"
                    ", labelList&, labelList&) const"
                )   << "in patch " << name()
                    << " : "
                    << "Cannot find point on face " << pp[oldFaceI]
                    << " with vertices "
                    << IndirectList<point>(pp.points(), pp[oldFaceI])()
                    << " that matches point " << wantedAnchor
                    << " when matching the halves of processor patch " << name()
                    << "Continuing with incorrect face ordering from now on!"
                    << endl;

                return false;
            }
        }

        ownerPatchPtr_.clear();

        // Return false if no change neccesary, true otherwise.

        forAll(faceMap, faceI)
        {
            if (faceMap[faceI] != faceI || rotation[faceI] != 0)
            {
                return true;
            }
        }

        return false;
    }
}


void Foam::cyclicPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    os.writeKeyword("neighbourPatch") << neighbPatchName_
        << token::END_STATEMENT << nl;
    switch (transform_)
    {
        case ROTATIONAL:
        {
            os.writeKeyword("transform") << transformTypeNames[ROTATIONAL]
                << token::END_STATEMENT << nl;
            os.writeKeyword("rotationAxis") << rotationAxis_
                << token::END_STATEMENT << nl;
            os.writeKeyword("rotationCentre") << rotationCentre_
                << token::END_STATEMENT << nl;
            break;
        }
        case TRANSLATIONAL:
        {
            os.writeKeyword("transform") << transformTypeNames[TRANSLATIONAL]
                << token::END_STATEMENT << nl;
            os.writeKeyword("separationVector") << separationVector_
                << token::END_STATEMENT << nl;
            break;
        }
        default:
        {
            // no additional info to write
        }
    }
}


// ************************************************************************* //
