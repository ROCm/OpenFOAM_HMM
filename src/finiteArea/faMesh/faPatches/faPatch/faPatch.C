/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "faPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "faBoundaryMesh.H"
#include "faMesh.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "polyMesh.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faPatch, 0);
    defineRunTimeSelectionTable(faPatch, dictionary);
    addToRunTimeSelectionTable(faPatch, faPatch, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faPatch::clearOut()
{
    deleteDemandDrivenData(edgeFacesPtr_);
    deleteDemandDrivenData(pointLabelsPtr_);
    deleteDemandDrivenData(pointEdgesPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faPatch::faPatch
(
    const word& name,
    const labelList& edgeLabels,
    const label index,
    const faBoundaryMesh& bm,
    const label ngbPolyPatchIndex
)
:
    labelList(edgeLabels),
    patchIdentifier(name, index),
    nbrPolyPatchId_(ngbPolyPatchIndex),
    boundaryMesh_(bm),
    edgeFacesPtr_(nullptr),
    pointLabelsPtr_(nullptr),
    pointEdgesPtr_(nullptr)
{}


Foam::faPatch::faPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const faBoundaryMesh& bm
)
:
    labelList(dict.get<labelList>("edgeLabels")),
    patchIdentifier(name, dict, index),
    nbrPolyPatchId_(dict.get<label>("ngbPolyPatchIndex")),
    boundaryMesh_(bm),
    edgeFacesPtr_(nullptr),
    pointLabelsPtr_(nullptr),
    pointEdgesPtr_(nullptr)
{}


Foam::faPatch::faPatch(const faPatch& p, const faBoundaryMesh& bm)
:
    labelList(p),
    patchIdentifier(p, p.index()),
    nbrPolyPatchId_(p.nbrPolyPatchId_),
    boundaryMesh_(bm),
    edgeFacesPtr_(nullptr),
    pointLabelsPtr_(nullptr),
    pointEdgesPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faPatch::~faPatch()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::faPatch::ngbPolyPatchIndex() const noexcept
{
    return nbrPolyPatchId_;
}


const Foam::faBoundaryMesh& Foam::faPatch::boundaryMesh() const noexcept
{
    return boundaryMesh_;
}


Foam::label Foam::faPatch::start() const
{
    return boundaryMesh().mesh().patchStarts()[index()];
}


const Foam::labelList& Foam::faPatch::pointLabels() const
{
    if (!pointLabelsPtr_)
    {
        calcPointLabels();
    }

    return *pointLabelsPtr_;
}


void Foam::faPatch::calcPointLabels() const
{
    SLList<label> labels;

    UList<edge> edges = patchSlice(boundaryMesh().mesh().edges());

    forAll(edges, edgeI)
    {
        bool existStart = false;
        bool existEnd = false;

        forAllIters(labels, iter)
        {
            if (*iter == edges[edgeI].start())
            {
                existStart = true;
            }

            if (*iter == edges[edgeI].end())
            {
                existEnd = true;
            }
        }

        if (!existStart)
        {
            labels.append(edges[edgeI].start());
        }

        if (!existEnd)
        {
            labels.append(edges[edgeI].end());
        }
    }

    pointLabelsPtr_ = new labelList(labels);
}


void Foam::faPatch::calcPointEdges() const
{
    const labelList& points = pointLabels();

    const edgeList::subList e = patchSlice(boundaryMesh().mesh().edges());

    // set up storage for pointEdges
    List<SLList<label>> pointEdgs(points.size());

    forAll(e, edgeI)
    {
        const edge& curPoints = e[edgeI];

        forAll(curPoints, pointI)
        {
            const label localPointIndex = points.find(curPoints[pointI]);

            pointEdgs[localPointIndex].append(edgeI);
        }
    }

    // sort out the list
    pointEdgesPtr_ = new labelListList(pointEdgs.size());
    labelListList& pEdges = *pointEdgesPtr_;

    forAll(pointEdgs, pointI)
    {
        pEdges[pointI].setSize(pointEdgs[pointI].size());

        label i = 0;
        for
        (
            SLList<label>::iterator curEdgesIter = pointEdgs[pointI].begin();
            curEdgesIter != pointEdgs[pointI].end();
            ++curEdgesIter, ++i
        )
        {
            pEdges[pointI][i] = curEdgesIter();
        }
    }
}


const Foam::labelListList& Foam::faPatch::pointEdges() const
{
    if (!pointEdgesPtr_)
    {
        calcPointEdges();
    }

    return *pointEdgesPtr_;
}


Foam::labelList Foam::faPatch::ngbPolyPatchFaces() const
{
    if (nbrPolyPatchId_ < 0)
    {
        return labelList();
    }

    labelList ngbFaces(faPatch::size());

    const faMesh& aMesh = boundaryMesh().mesh();
    const polyMesh& pMesh = aMesh.mesh();
    const auto& patch = aMesh.patch();

    const labelListList& edgeFaces = pMesh.edgeFaces();

    const labelList meshEdges
    (
        patch.meshEdges(pMesh.edges(), pMesh.pointEdges())
    );

    forAll(ngbFaces, edgeI)
    {
        ngbFaces[edgeI] = -1;

        label curEdge = (*this)[edgeI];

        label curPMeshEdge = meshEdges[curEdge];

        forAll(edgeFaces[curPMeshEdge], faceI)
        {
            label curFace = edgeFaces[curPMeshEdge][faceI];

            label curPatchID = pMesh.boundaryMesh().whichPatch(curFace);

            if (curPatchID == nbrPolyPatchId_)
            {
                ngbFaces[edgeI] = curFace;
            }
        }

        if (ngbFaces[edgeI] == -1)
        {
            WarningInFunction
                << "Problem with determination of edge ngb faces!"
                << endl;
        }
    }

    return ngbFaces;
}


Foam::tmp<Foam::vectorField> Foam::faPatch::ngbPolyPatchFaceNormals() const
{
    if (nbrPolyPatchId_ < 0)
    {
        return tmp<vectorField>::New();
    }

    auto tfN = tmp<vectorField>::New();
    auto& fN = tfN.ref();

    fN.setSize(faPatch::size());

    labelList ngbFaces = ngbPolyPatchFaces();

    const polyMesh& pMesh = boundaryMesh().mesh()();

    const faceList& faces = pMesh.faces();
    const pointField& points = pMesh.points();

    forAll(fN, faceI)
    {
        fN[faceI] = faces[ngbFaces[faceI]].unitNormal(points);
    }

    return tfN;
}


Foam::tmp<Foam::vectorField> Foam::faPatch::ngbPolyPatchPointNormals() const
{
    if (nbrPolyPatchId_ < 0)
    {
        return tmp<vectorField>::New();
    }

    const labelListList& pntEdges = pointEdges();

    auto tpN = tmp<vectorField>::New(pntEdges.size(), Zero);
    auto& pN = tpN.ref();

    const vectorField faceNormals(ngbPolyPatchFaceNormals());

    forAll(pN, pointI)
    {
        forAll(pntEdges[pointI], edgeI)
        {
            pN[pointI] += faceNormals[pntEdges[pointI][edgeI]];
        }
    }

    pN /= mag(pN);

    return tpN;
}


const Foam::labelUList& Foam::faPatch::edgeFaces() const
{
    if (!edgeFacesPtr_)
    {
        edgeFacesPtr_ = new labelList::subList
        (
            patchSlice(boundaryMesh().mesh().edgeOwner())
        );
    }

    return *edgeFacesPtr_;
}


const Foam::vectorField& Foam::faPatch::edgeCentres() const
{
    return boundaryMesh().mesh().edgeCentres().boundaryField()[index()];
}


const Foam::vectorField& Foam::faPatch::edgeLengths() const
{
    return boundaryMesh().mesh().Le().boundaryField()[index()];
}


const Foam::scalarField& Foam::faPatch::magEdgeLengths() const
{
    return boundaryMesh().mesh().magLe().boundaryField()[index()];
}


Foam::tmp<Foam::vectorField> Foam::faPatch::edgeNormals() const
{
    tmp<vectorField> eN(new vectorField(size()));

    eN.ref() = edgeLengths()/magEdgeLengths();

    return eN;
}


Foam::tmp<Foam::vectorField> Foam::faPatch::edgeFaceCentres() const
{
    auto tfc = tmp<vectorField>::New(size());
    auto& fc = tfc.ref();

    // get reference to global face centres
    const vectorField& gfc =
        boundaryMesh().mesh().areaCentres().internalField();

    const labelUList& faceLabels = edgeFaces();

    forAll(faceLabels, edgeI)
    {
        fc[edgeI] = gfc[faceLabels[edgeI]];
    }

    return tfc;
}


Foam::tmp<Foam::vectorField> Foam::faPatch::delta() const
{
    return edgeNormals()*(edgeNormals() & (edgeCentres() - edgeFaceCentres()));
}


void Foam::faPatch::makeDeltaCoeffs(scalarField& dc) const
{
    dc = scalar(1)/(edgeNormals() & delta());
}


void Foam::faPatch::makeCorrectionVectors(vectorField& k) const
{
    vectorField unitDelta(delta()/mag(delta()));
    vectorField edgeNormMag(edgeNormals()/mag(edgeNormals()));
    scalarField dn(edgeNormals() & delta());

    k = edgeNormMag - (scalar(1)/(unitDelta & edgeNormMag))*unitDelta;
}


const Foam::scalarField& Foam::faPatch::deltaCoeffs() const
{
    return boundaryMesh().mesh().deltaCoeffs().boundaryField()[index()];
}


void Foam::faPatch::makeWeights(scalarField& w) const
{
    w = scalar(1);
}


const Foam::scalarField& Foam::faPatch::weights() const
{
    return boundaryMesh().mesh().weights().boundaryField()[index()];
}


void Foam::faPatch::movePoints(const pointField& points)
{}


void Foam::faPatch::resetEdges(const UList<label>& newEdges)
{
    clearOut();
    static_cast<labelList&>(*this) = newEdges;
}


void Foam::faPatch::resetEdges(labelList&& newEdges)
{
    clearOut();
    static_cast<labelList&>(*this) = std::move(newEdges);
}


void Foam::faPatch::write(Ostream& os) const
{
    os.writeEntry("type", type());

    patchIdentifier::write(os);

    os.writeEntry("ngbPolyPatchIndex", nbrPolyPatchId_);
    static_cast<const labelList&>(*this).writeEntry("edgeLabels", os);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const faPatch& p)
{
    p.write(os);
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
