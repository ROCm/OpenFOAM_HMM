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
#include "edgeHashes.H"
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
    patchIdentifier(name, index),
    labelList(edgeLabels),
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
    patchIdentifier(name, dict, index),
    labelList(dict.get<labelList>("edgeLabels")),
    nbrPolyPatchId_(dict.get<label>("ngbPolyPatchIndex")),
    boundaryMesh_(bm),
    edgeFacesPtr_(nullptr),
    pointLabelsPtr_(nullptr),
    pointEdgesPtr_(nullptr)
{}


Foam::faPatch::faPatch(const faPatch& p, const faBoundaryMesh& bm)
:
    patchIdentifier(p, p.index()),
    labelList(p),
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

const Foam::faBoundaryMesh& Foam::faPatch::boundaryMesh() const noexcept
{
    return boundaryMesh_;
}


Foam::label Foam::faPatch::start() const
{
    return boundaryMesh().mesh().patchStarts()[index()];
}


Foam::List<Foam::labelPair> Foam::faPatch::boundaryConnections() const
{
    const auto& connections = boundaryMesh().mesh().boundaryConnections();
    const label nInternalEdges = boundaryMesh().mesh().nInternalEdges();

    List<labelPair> output(this->nEdges());

    // Like an IndirectList but removing the nInternalEdges offset
    label count = 0;
    for (const label patchEdgei : this->edgeLabels())
    {
        const label bndEdgei = (patchEdgei - nInternalEdges);
        output[count] = connections[bndEdgei];
        ++count;
    }

    return output;
}


Foam::labelList Foam::faPatch::boundaryProcs() const
{
    const auto& connections = boundaryMesh().mesh().boundaryConnections();
    const label nInternalEdges = boundaryMesh().mesh().nInternalEdges();

    labelHashSet procsUsed(2*Pstream::nProcs());

    for (const label patchEdgei : this->edgeLabels())
    {
        const label bndEdgei = (patchEdgei - nInternalEdges);
        const label proci = connections[bndEdgei].first();
        procsUsed.insert(proci);
    }
    procsUsed.erase(-1);  // placeholder value
    procsUsed.erase(Pstream::myProcNo());

    return procsUsed.sortedToc();
}


Foam::List<Foam::labelPair> Foam::faPatch::boundaryProcSizes() const
{
    const auto& connections = boundaryMesh().mesh().boundaryConnections();
    const label nInternalEdges = boundaryMesh().mesh().nInternalEdges();

    Map<label> procCount(2*Pstream::nProcs());

    for (const label patchEdgei : this->edgeLabels())
    {
        const label bndEdgei = (patchEdgei - nInternalEdges);
        const label proci = connections[bndEdgei].first();
        ++procCount(proci);
    }
    procCount.erase(-1);  // placeholder value
    procCount.erase(Pstream::myProcNo());

    // Flatten as list
    List<labelPair> output(procCount.size());
    label count = 0;
    for (const label proci : procCount.sortedToc())
    {
        output[count].first() = proci;
        output[count].second() = procCount[proci];  // size
        ++count;
    }

    return output;
}


const Foam::labelList& Foam::faPatch::pointLabels() const
{
    if (!pointLabelsPtr_)
    {
        calcPointLabels();
    }

    return *pointLabelsPtr_;
}


const Foam::labelListList& Foam::faPatch::pointEdges() const
{
    if (!pointEdgesPtr_)
    {
        calcPointEdges();
    }

    return *pointEdgesPtr_;
}


void Foam::faPatch::calcPointLabels() const
{
    const edgeList::subList edges = patchSlice(boundaryMesh().mesh().edges());

    // Walk boundary edges.
    // The edge orientation corresponds to the face orientation
    // (outwards normal).

    // Note: could combine this with calcPointEdges for more efficiency

    // Map<label> markedPoints(4*edges.size());
    labelHashSet markedPoints(4*edges.size());
    DynamicList<label> dynEdgePoints(2*edges.size());

    for (const edge& e : edges)
    {
        // if (markedPoints.insert(e.first(), markedPoints.size()))
        if (markedPoints.insert(e.first()))
        {
            dynEdgePoints.append(e.first());
        }
        // if (markedPoints.insert(e.second(), markedPoints.size()))
        if (markedPoints.insert(e.second()))
        {
            dynEdgePoints.append(e.second());
        }
    }

    // Transfer to plain list (reuse storage)
    pointLabelsPtr_ = new labelList(std::move(dynEdgePoints));
    /// const labelList& edgePoints = *pointLabelsPtr_;
    ///
    /// // Cannot use invertManyToMany - we have non-local edge numbering
    ///
    /// // Intermediate storage for pointEdges.
    /// // Points on the boundary will normally connect 1 or 2 edges only.
    /// List<DynamicList<label,2>> dynPointEdges(edgePoints.size());
    ///
    /// forAll(edges, edgei)
    /// {
    ///     const edge& e = edges[edgei];
    ///
    ///     dynPointEdges[markedPoints[e.first()]].append(edgei);
    ///     dynPointEdges[markedPoints[e.second()]].append(edgei);
    /// }
    ///
    /// // Flatten to regular list
    /// pointEdgesPtr_ = new labelListList(edgePoints.size());
    /// auto& pEdges = *pointEdgesPtr_;
    ///
    /// forAll(pEdges, pointi)
    /// {
    ///     pEdges[pointi] = std::move(dynPointEdges[pointi]);
    /// }
}


void Foam::faPatch::calcPointEdges() const
{
    const edgeList::subList edges = patchSlice(boundaryMesh().mesh().edges());

    const labelList& edgePoints = pointLabels();

    // Cannot use invertManyToMany - we have non-local edge numbering

    // Intermediate storage for pointEdges.
    // Points on the boundary will normally connect 1 or 2 edges only.
    List<DynamicList<label,2>> dynPointEdges(edgePoints.size());

    forAll(edges, edgei)
    {
        const edge& e = edges[edgei];

        dynPointEdges[edgePoints.find(e.first())].append(edgei);
        dynPointEdges[edgePoints.find(e.second())].append(edgei);
    }

    // Flatten to regular list
    pointEdgesPtr_ = new labelListList(edgePoints.size());
    auto& pEdges = *pointEdgesPtr_;

    forAll(pEdges, pointi)
    {
        pEdges[pointi] = std::move(dynPointEdges[pointi]);
    }
}


Foam::tmp<Foam::vectorField> Foam::faPatch::ngbPolyPatchFaceNormals() const
{
    if (nbrPolyPatchId_ < 0)
    {
        return tmp<vectorField>::New();
    }

    return boundaryMesh().mesh().haloFaceNormals(this->index());
}


Foam::tmp<Foam::vectorField> Foam::faPatch::ngbPolyPatchPointNormals() const
{
    if (nbrPolyPatchId_ < 0)
    {
        return tmp<vectorField>::New();
    }

    // Unit normals for the neighbour patch faces
    const vectorField faceNormals
    (
        boundaryMesh().mesh().haloFaceNormals(this->index())
    );

    const labelListList& pntEdges = pointEdges();

    auto tpointNorm = tmp<vectorField>::New(pntEdges.size());
    auto& pointNorm = tpointNorm.ref();

    forAll(pointNorm, pointi)
    {
        vector& n = pointNorm[pointi];
        n = Zero;

        for (const label bndEdgei : pntEdges[pointi])
        {
            n += faceNormals[bndEdgei];
        }

        n.normalise();
    }

    return tpointNorm;
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
    auto tedgeNorm = tmp<vectorField>::New(edgeLengths());

    for (vector& n : tedgeNorm.ref())
    {
        n.normalise();
    }

    return tedgeNorm;
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
