/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "turbulentDFSEMInletFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "momentOfInertia.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::label Foam::turbulentDFSEMInletFvPatchVectorField::seedIterMax_ = 1000;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::turbulentDFSEMInletFvPatchVectorField::writeEddyOBJ() const
{
    {
        // Output the bounding box
        OFstream os(db().time().path()/"eddyBox.obj");

        const polyPatch& pp = this->patch().patch();
        const labelList& boundaryPoints = pp.boundaryPoints();
        const pointField& localPoints = pp.localPoints();

        const vector offset(patchNormal_*maxSigmaX_);
        forAll(boundaryPoints, i)
        {
            point p = localPoints[boundaryPoints[i]];
            p += offset;
            os  << "v " << p.x() << " " << p.y() << " " << p.z() << nl;
        }

        forAll(boundaryPoints, i)
        {
            point p = localPoints[boundaryPoints[i]];
            p -= offset;
            os  << "v " << p.x() << " " << p.y() << " " << p.z() << nl;
        }
    }

    {
        const Time& time = db().time();
        OFstream os
        (
            time.path()/"eddies_" + Foam::name(time.timeIndex()) + ".obj"
        );

        label pointOffset = 0;
        forAll(eddies_, eddyI)
        {
            const eddy& e = eddies_[eddyI];
            pointOffset += e.writeSurfaceOBJ(pointOffset, patchNormal_, os);
        }
    }
}


void Foam::turbulentDFSEMInletFvPatchVectorField::writeLumleyCoeffs() const
{
    // Output list of xi vs eta

    OFstream os(db().time().path()/"lumley_interpolated.out");

    os  << "# xi" << token::TAB << "eta" << endl;

    const scalar t = db().time().timeOutputValue();
    const symmTensorField R(R_->value(t)/sqr(Uref_));

    forAll(R, faceI)
    {
        // Normalised anisotropy tensor
        const symmTensor devR(dev(R[faceI]/(tr(R[faceI]))));

        // Second tensor invariant
        const scalar ii = min(0, invariantII(devR));

        // Third tensor invariant
        const scalar iii = invariantIII(devR);

        // xi, eta
        // See Pope - characterization of Reynolds-stress anisotropy
        const scalar xi = cbrt(0.5*iii);
        const scalar eta = sqrt(-ii/3.0);
        os  << xi << token::TAB << eta << token::TAB
            << ii << token::TAB << iii << endl;
    }
}


void Foam::turbulentDFSEMInletFvPatchVectorField::initialisePatch()
{
    const vectorField nf(patch().nf());

    // Patch normal points into domain
    patchNormal_ = -gAverage(nf);

    // Check that patch is planar
    const scalar error = max(magSqr(patchNormal_ + nf));

    if (error > SMALL)
    {
        WarningInFunction
            << "Patch " << patch().name() << " is not planar"
            << endl;
    }

    patchNormal_ /= mag(patchNormal_) + ROOTVSMALL;


    // Decompose the patch faces into triangles to enable point search

    const polyPatch& patch = this->patch().patch();
    const pointField& points = patch.points();

    // Triangulate the patch faces and create addressing
    DynamicList<label> triToFace(2*patch.size());
    DynamicList<scalar> triMagSf(2*patch.size());
    DynamicList<face> triFace(2*patch.size());
    DynamicList<face> tris(5);

    // Set zero value at the start of the tri area list
    triMagSf.append(0.0);

    forAll(patch, faceI)
    {
        const face& f = patch[faceI];

        tris.clear();
        f.triangles(points, tris);

        forAll(tris, i)
        {
            triToFace.append(faceI);
            triFace.append(tris[i]);
            triMagSf.append(tris[i].mag(points));
        }
    }

    for (auto& s : sumTriMagSf_)
    {
        s = 0.0;
    }

    sumTriMagSf_[Pstream::myProcNo() + 1] = sum(triMagSf);

    Pstream::listCombineGather(sumTriMagSf_, maxEqOp<scalar>());
    Pstream::listCombineScatter(sumTriMagSf_);

    for (label i = 1; i < triMagSf.size(); ++i)
    {
        triMagSf[i] += triMagSf[i-1];
    }

    // Transfer to persistent storage
    triFace_.transfer(triFace);
    triToFace_.transfer(triToFace);
    triCumulativeMagSf_.transfer(triMagSf);

    // Convert sumTriMagSf_ into cumulative sum of areas per proc
    for (label i = 1; i < sumTriMagSf_.size(); ++i)
    {
        sumTriMagSf_[i] += sumTriMagSf_[i-1];
    }

    // Global patch area (over all processors)
    patchArea_ = sumTriMagSf_.last();

    // Local patch bounds (this processor)
    patchBounds_ = boundBox(patch.localPoints(), false);
    patchBounds_.inflate(0.1);

    // Determine if all eddies spawned from a single processor
    singleProc_ = patch.size() == returnReduce(patch.size(), sumOp<label>());
    reduce(singleProc_, orOp<bool>());
}


void Foam::turbulentDFSEMInletFvPatchVectorField::initialiseEddyBox()
{
    const scalarField& magSf = patch().magSf();

    const scalarField L(L_->value(db().time().timeOutputValue())/Lref_);

    // (PCF:Eq. 14)
    const scalarField cellDx(max(Foam::sqrt(magSf), 2/patch().deltaCoeffs()));

    // Inialise eddy box extents
    forAll(*this, faceI)
    {
        scalar& s = sigmax_[faceI];

        // Average length scale (SST:Eq. 24)
        // Personal communication regarding (PCR:Eq. 14)
        //  - the min operator in Eq. 14 is a typo, and should be a max operator
        s = min(mag(L[faceI]), kappa_*delta_);
        s = max(s, nCellPerEddy_*cellDx[faceI]);
    }

    // Maximum extent across all processors
    maxSigmaX_ = gMax(sigmax_);

    // Eddy box volume
    v0_ = 2*gSum(magSf)*maxSigmaX_;

    if (debug)
    {
        Info<< "Patch: " << patch().patch().name() << " eddy box:" << nl
            << "    volume    : " << v0_ << nl
            << "    maxSigmaX : " << maxSigmaX_ << nl
            << endl;
    }
}


Foam::pointIndexHit Foam::turbulentDFSEMInletFvPatchVectorField::setNewPosition
(
    const bool global
)
{
    // Initialise to miss
    pointIndexHit pos(false, vector::max, -1);

    const polyPatch& patch = this->patch().patch();
    const pointField& points = patch.points();

    if (global)
    {
        const scalar areaFraction =
            rndGen_.globalPosition<scalar>(0, patchArea_);

        // Determine which processor to use
        label procI = 0;
        forAllReverse(sumTriMagSf_, i)
        {
            if (areaFraction >= sumTriMagSf_[i])
            {
                procI = i;
                break;
            }
        }

        if (Pstream::myProcNo() == procI)
        {
            // Find corresponding decomposed face triangle
            label triI = 0;
            const scalar offset = sumTriMagSf_[procI];
            forAllReverse(triCumulativeMagSf_, i)
            {
                if (areaFraction > triCumulativeMagSf_[i] + offset)
                {
                    triI = i;
                    break;
                }
            }

            // Find random point in triangle
            const face& tf = triFace_[triI];
            const triPointRef tri(points[tf[0]], points[tf[1]], points[tf[2]]);

            pos.setHit();
            pos.setIndex(triToFace_[triI]);
            pos.rawPoint() = tri.randomPoint(rndGen_);
        }
    }
    else
    {
        // Find corresponding decomposed face triangle on local processor
        label triI = 0;
        const scalar maxAreaLimit = triCumulativeMagSf_.last();
        const scalar areaFraction = rndGen_.position<scalar>(0, maxAreaLimit);

        forAllReverse(triCumulativeMagSf_, i)
        {
            if (areaFraction > triCumulativeMagSf_[i])
            {
                triI = i;
                break;
            }
        }

        // Find random point in triangle
        const face& tf = triFace_[triI];
        const triPointRef tri(points[tf[0]], points[tf[1]], points[tf[2]]);

        pos.setHit();
        pos.setIndex(triToFace_[triI]);
        pos.rawPoint() = tri.randomPoint(rndGen_);
    }

    return pos;
}


void Foam::turbulentDFSEMInletFvPatchVectorField::initialiseEddies()
{
    const scalar t = db().time().timeOutputValue();
    const symmTensorField R(R_->value(t)/sqr(Uref_));

    DynamicList<eddy> eddies(size());

    // Initialise eddy properties
    scalar sumVolEddy = 0;
    scalar sumVolEddyAllProc = 0;

    while (sumVolEddyAllProc/v0_ < d_)
    {
        bool search = true;
        label iter = 0;

        while (search && iter++ < seedIterMax_)
        {
            // Get new parallel consistent position
            pointIndexHit pos(setNewPosition(true));
            const label patchFaceI = pos.index();

            // Note: only 1 processor will pick up this face
            if (patchFaceI != -1)
            {
                eddy e
                (
                    patchFaceI,
                    pos.hitPoint(),
                    rndGen_.position<scalar>(-maxSigmaX_, maxSigmaX_),
                    sigmax_[patchFaceI],
                    R[patchFaceI],
                    rndGen_
                );

                // If eddy valid, patchFaceI is non-zero
                if (e.patchFaceI() != -1)
                {
                    eddies.append(e);
                    sumVolEddy += e.volume();
                    search = false;
                }
            }
            // else eddy on remote processor

            reduce(search, andOp<bool>());
        }


        sumVolEddyAllProc = returnReduce(sumVolEddy, sumOp<scalar>());
    }
    eddies_.transfer(eddies);

    nEddy_ = eddies_.size();

    if (debug)
    {
        Pout<< "Patch:" << patch().patch().name();

        if (Pstream::parRun())
        {
            Pout<< " processor:" << Pstream::myProcNo();
        }

        Pout<< " seeded:" << nEddy_ << " eddies" << endl;
    }

    reduce(nEddy_, sumOp<label>());

    if (nEddy_ > 0)
    {
        Info<< "Turbulent DFSEM patch: " << patch().name()
            << " seeded " << nEddy_ << " eddies with total volume "
            << sumVolEddyAllProc
            << endl;
    }
    else
    {
        WarningInFunction
            << "Patch: " << patch().patch().name()
            << " on field " << internalField().name()
            << ": No eddies seeded - please check your set-up"
            << endl;
    }
}


void Foam::turbulentDFSEMInletFvPatchVectorField::convectEddies
(
    const vector& UBulk,
    const scalar deltaT
)
{
    const scalar t = db().time().timeOutputValue();
    const symmTensorField R(R_->value(t)/sqr(Uref_));

    // Note: all operations applied to local processor only

    label nRecycled = 0;

    forAll(eddies_, eddyI)
    {
        eddy& e = eddies_[eddyI];
        e.move(deltaT*(UBulk & patchNormal_));

        const scalar position0 = e.x();

        // Check to see if eddy has exited downstream box plane
        if (position0 > maxSigmaX_)
        {
            bool search = true;
            label iter = 0;

            while (search && iter++ < seedIterMax_)
            {
                // Spawn new eddy with new random properties (intensity etc)
                pointIndexHit pos(setNewPosition(false));
                const label patchFaceI = pos.index();

                e = eddy
                    (
                        patchFaceI,
                        pos.hitPoint(),
                        position0 - floor(position0/maxSigmaX_)*maxSigmaX_,
                        sigmax_[patchFaceI],
                        R[patchFaceI],
                        rndGen_
                    );

                if (e.patchFaceI() != -1)
                {
                    search = false;
                }
            }

            nRecycled++;
        }
    }

    reduce(nRecycled, sumOp<label>());

    if (debug && nRecycled > 0)
    {
        Info<< "Patch: " << patch().patch().name()
            << " recycled " << nRecycled << " eddies"
            << endl;
    }
}


Foam::vector Foam::turbulentDFSEMInletFvPatchVectorField::uPrimeEddy
(
    const List<eddy>& eddies,
    const point& patchFaceCf
) const
{
    vector uPrime(Zero);

    forAll(eddies, k)
    {
        const eddy& e = eddies[k];
        uPrime += e.uPrime(patchFaceCf, patchNormal_);
    }

    return uPrime;
}


void Foam::turbulentDFSEMInletFvPatchVectorField::calcOverlappingProcEddies
(
    List<List<eddy>>& overlappingEddies
) const
{
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    List<boundBox> patchBBs(Pstream::nProcs());
    patchBBs[Pstream::myProcNo()] = patchBounds_;
    Pstream::gatherList(patchBBs);
    Pstream::scatterList(patchBBs);

    // Per processor indices into all segments to send
    List<DynamicList<label>> dynSendMap(Pstream::nProcs());

    forAll(eddies_, i)
    {
        // Collect overlapping eddies
        const eddy& e = eddies_[i];

        // Eddy bounds
        const point x(e.position(patchNormal_));
        boundBox ebb(e.bounds());
        ebb.min() += x;
        ebb.max() += x;

        forAll(patchBBs, procI)
        {
            // Not including intersection with local patch
            if (procI != Pstream::myProcNo())
            {
                if (ebb.overlaps(patchBBs[procI]))
                {
                    dynSendMap[procI].append(i);
                }
            }
        }
    }

    labelListList sendMap(Pstream::nProcs());
    forAll(sendMap, procI)
    {
        sendMap[procI].transfer(dynSendMap[procI]);
    }

    // Send the number of eddies for local processors to receive
    labelListList sendSizes(Pstream::nProcs());
    sendSizes[Pstream::myProcNo()].setSize(Pstream::nProcs());
    forAll(sendMap, procI)
    {
        sendSizes[Pstream::myProcNo()][procI] = sendMap[procI].size();
    }
    Pstream::gatherList(sendSizes);
    Pstream::scatterList(sendSizes);

    // Determine order of receiving
    labelListList constructMap(Pstream::nProcs());

    // Local segment first
    constructMap[Pstream::myProcNo()] = identity
    (
        sendMap[Pstream::myProcNo()].size()
    );

    label segmentI = constructMap[Pstream::myProcNo()].size();
    forAll(constructMap, procI)
    {
        if (procI != Pstream::myProcNo())
        {
            // What I need to receive is what other processor is sending to me
            const label nRecv = sendSizes[procI][Pstream::myProcNo()];
            constructMap[procI].setSize(nRecv);

            for (label i = 0; i < nRecv; ++i)
            {
                constructMap[procI][i] = segmentI++;
            }
        }
    }

    mapDistribute map(segmentI, std::move(sendMap), std::move(constructMap));

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    for (const int domain : Pstream::allProcs())
    {
        const labelList& sendElems = map.subMap()[domain];

        if (domain != Pstream::myProcNo() && sendElems.size())
        {
            List<eddy> subEddies(UIndirectList<eddy>(eddies_, sendElems));

            UOPstream toDomain(domain, pBufs);

            toDomain<< subEddies;
        }
    }

    // Start receiving
    pBufs.finishedSends();

    // Consume
    for (const int domain : Pstream::allProcs())
    {
        const labelList& recvElems = map.constructMap()[domain];

        if (domain != Pstream::myProcNo() && recvElems.size())
        {
            UIPstream str(domain, pBufs);
            {
                str >> overlappingEddies[domain];
            }
        }
    }

    // Restore tag
    UPstream::msgType() = oldTag;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentDFSEMInletFvPatchVectorField::
turbulentDFSEMInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    U_(nullptr),
    R_(nullptr),
    L_(nullptr),
    delta_(1.0),
    d_(1.0),
    kappa_(0.41),
    Uref_(1.0),
    Lref_(1.0),
    scale_(1.0),
    m_(0.5),
    nCellPerEddy_(5),

    patchArea_(-1),
    triFace_(),
    triToFace_(),
    triCumulativeMagSf_(),
    sumTriMagSf_(Pstream::nProcs() + 1, Zero),
    patchNormal_(Zero),
    patchBounds_(boundBox::invertedBox),

    eddies_(Zero),
    v0_(Zero),
    rndGen_(Pstream::myProcNo()),
    sigmax_(size(), Zero),
    maxSigmaX_(Zero),
    nEddy_(0),
    curTimeIndex_(-1),
    singleProc_(false),
    writeEddies_(false)
{}


Foam::turbulentDFSEMInletFvPatchVectorField::
turbulentDFSEMInletFvPatchVectorField
(
    const turbulentDFSEMInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    U_(ptf.U_.clone(patch().patch())),
    R_(ptf.R_.clone(patch().patch())),
    L_(ptf.L_.clone(patch().patch())),
    delta_(ptf.delta_),
    d_(ptf.d_),
    kappa_(ptf.kappa_),
    Uref_(ptf.Uref_),
    Lref_(ptf.Lref_),
    scale_(ptf.scale_),
    m_(ptf.m_),
    nCellPerEddy_(ptf.nCellPerEddy_),

    patchArea_(ptf.patchArea_),
    triFace_(ptf.triFace_),
    triToFace_(ptf.triToFace_),
    triCumulativeMagSf_(ptf.triCumulativeMagSf_),
    sumTriMagSf_(ptf.sumTriMagSf_),
    patchNormal_(ptf.patchNormal_),
    patchBounds_(ptf.patchBounds_),

    eddies_(ptf.eddies_),
    v0_(ptf.v0_),
    rndGen_(ptf.rndGen_),
    sigmax_(ptf.sigmax_, mapper),
    maxSigmaX_(ptf.maxSigmaX_),
    nEddy_(ptf.nEddy_),
    curTimeIndex_(-1),
    singleProc_(ptf.singleProc_),
    writeEddies_(ptf.writeEddies_)
{}


Foam::turbulentDFSEMInletFvPatchVectorField::
turbulentDFSEMInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    U_(PatchFunction1<vector>::New(patch().patch(), "U", dict)),
    R_(PatchFunction1<symmTensor>::New(patch().patch(), "R", dict)),
    L_(PatchFunction1<scalar>::New(patch().patch(), "L", dict)),
    delta_(dict.getCheck<scalar>("delta", scalarMinMax::ge(0))),
    d_(dict.getCheckOrDefault<scalar>("d", 1, scalarMinMax::ge(SMALL))),
    kappa_(dict.getCheckOrDefault<scalar>("kappa", 0.41, scalarMinMax::ge(0))),
    Uref_(dict.getCheckOrDefault<scalar>("Uref", 1, scalarMinMax::ge(SMALL))),
    Lref_(dict.getCheckOrDefault<scalar>("Lref", 1, scalarMinMax::ge(SMALL))),
    scale_(dict.getCheckOrDefault<scalar>("scale", 1, scalarMinMax::ge(0))),
    m_(dict.getCheckOrDefault<scalar>("m", 0.5, scalarMinMax::ge(0))),
    nCellPerEddy_(dict.getOrDefault<label>("nCellPerEddy", 5)),

    patchArea_(-1),
    triFace_(),
    triToFace_(),
    triCumulativeMagSf_(),
    sumTriMagSf_(Pstream::nProcs() + 1, Zero),
    patchNormal_(Zero),
    patchBounds_(boundBox::invertedBox),

    eddies_(),
    v0_(Zero),
    rndGen_(),
    sigmax_(size(), Zero),
    maxSigmaX_(Zero),
    nEddy_(0),
    curTimeIndex_(-1),
    singleProc_(false),
    writeEddies_(dict.getOrDefault("writeEddies", false))
{
    eddy::debug = debug;

    const scalar t = db().time().timeOutputValue();
    const symmTensorField R(R_->value(t)/sqr(Uref_));

    checkStresses(R);
}


Foam::turbulentDFSEMInletFvPatchVectorField::
turbulentDFSEMInletFvPatchVectorField
(
    const turbulentDFSEMInletFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    U_(ptf.U_.clone(patch().patch())),
    R_(ptf.R_.clone(patch().patch())),
    L_(ptf.L_.clone(patch().patch())),
    delta_(ptf.delta_),
    d_(ptf.d_),
    kappa_(ptf.kappa_),
    Uref_(ptf.Uref_),
    Lref_(ptf.Lref_),
    scale_(ptf.scale_),
    m_(ptf.m_),
    nCellPerEddy_(ptf.nCellPerEddy_),

    patchArea_(ptf.patchArea_),
    triFace_(ptf.triFace_),
    triToFace_(ptf.triToFace_),
    triCumulativeMagSf_(ptf.triCumulativeMagSf_),
    sumTriMagSf_(ptf.sumTriMagSf_),
    patchNormal_(ptf.patchNormal_),
    patchBounds_(ptf.patchBounds_),

    eddies_(ptf.eddies_),
    v0_(ptf.v0_),
    rndGen_(ptf.rndGen_),
    sigmax_(ptf.sigmax_),
    maxSigmaX_(ptf.maxSigmaX_),
    nEddy_(ptf.nEddy_),
    curTimeIndex_(-1),
    singleProc_(ptf.singleProc_),
    writeEddies_(ptf.writeEddies_)
{}


Foam::turbulentDFSEMInletFvPatchVectorField::
turbulentDFSEMInletFvPatchVectorField
(
    const turbulentDFSEMInletFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    U_(ptf.U_.clone(patch().patch())),
    R_(ptf.R_.clone(patch().patch())),
    L_(ptf.L_.clone(patch().patch())),
    delta_(ptf.delta_),
    d_(ptf.d_),
    kappa_(ptf.kappa_),
    Uref_(ptf.Uref_),
    Lref_(ptf.Lref_),
    scale_(ptf.scale_),
    m_(ptf.m_),
    nCellPerEddy_(ptf.nCellPerEddy_),

    patchArea_(ptf.patchArea_),
    triFace_(ptf.triFace_),
    triToFace_(ptf.triToFace_),
    triCumulativeMagSf_(ptf.triCumulativeMagSf_),
    sumTriMagSf_(ptf.sumTriMagSf_),
    patchNormal_(ptf.patchNormal_),
    patchBounds_(ptf.patchBounds_),

    eddies_(ptf.eddies_),
    v0_(ptf.v0_),
    rndGen_(ptf.rndGen_),
    sigmax_(ptf.sigmax_),
    maxSigmaX_(ptf.maxSigmaX_),
    nEddy_(ptf.nEddy_),
    curTimeIndex_(-1),
    singleProc_(ptf.singleProc_),
    writeEddies_(ptf.writeEddies_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::turbulentDFSEMInletFvPatchVectorField::checkStresses
(
    const symmTensorField& Rf
)
{
    // Perform checks of the stress tensor based on Cholesky decomposition
    // constraints

    forAll(Rf, facei)
    {
        const symmTensor& R = Rf[facei];

        if (R.xx() <= 0)
        {
            FatalErrorInFunction
                << "Reynolds stress " << R << " at face " << facei
                << " does not obey the constraint: R_xx > 0"
                << exit(FatalError);
        }

        const scalar a_xx = sqrt(R.xx());

        const scalar a_xy = R.xy()/a_xx;

        const scalar a_yy_2 = R.yy() - sqr(a_xy);

        if (a_yy_2 < 0)
        {
            FatalErrorInFunction
                << "Reynolds stress " << R << " at face " << facei
                << " leads to an invalid Cholesky decomposition due to the "
                << "constraint R_yy - sqr(a_xy) >= 0"
                << exit(FatalError);
        }

        const scalar a_yy = Foam::sqrt(a_yy_2);

        const scalar a_xz = R.xz()/a_xx;

        const scalar a_yz = (R.yz() - a_xy*a_xz)/a_yy;

        const scalar a_zz_2 = R.zz() - sqr(a_xz) - sqr(a_yz);

        if (a_zz_2 < 0)
        {
            FatalErrorInFunction
                << "Reynolds stress " << R << " at face " << facei
                << " leads to an invalid Cholesky decomposition due to the "
                << "constraint R_zz - sqr(a_xz) - sqr(a_yz) >= 0"
                << exit(FatalError);
        }

        const scalar a_zz = Foam::sqrt(a_zz_2);

        if (debug)
        {
            Pout<< "R: " << R
                << " a_xx:" << a_xx << " a_xy:" << a_xy << " a_xz:" << a_xy
                <<                     " a_yy:" << a_yy << " a_yz:" << a_yz
                <<                                         " a_zz:" << a_zz
                << endl;
        }
    }

    return true;
}


void Foam::turbulentDFSEMInletFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<vector>::autoMap(m);

    if (U_)
    {
        U_->autoMap(m);
    }
    if (R_)
    {
        R_->autoMap(m);
    }
    if (L_)
    {
        L_->autoMap(m);
    }

    sigmax_.autoMap(m);
}


void Foam::turbulentDFSEMInletFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<vector>::rmap(ptf, addr);

    const auto& dfsemptf =
        refCast<const turbulentDFSEMInletFvPatchVectorField>(ptf);

    if (U_)
    {
        U_->rmap(dfsemptf.U_(), addr);
    }
    if (R_)
    {
        R_->rmap(dfsemptf.R_(), addr);
    }
    if (L_)
    {
        L_->rmap(dfsemptf.L_(), addr);
    }

    sigmax_.rmap(dfsemptf.sigmax_, addr);
}


void Foam::turbulentDFSEMInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (curTimeIndex_ == -1)
    {
        initialisePatch();

        initialiseEddyBox();

        initialiseEddies();
    }


    if (curTimeIndex_ != db().time().timeIndex())
    {
        tmp<vectorField> UMean =
            U_->value(db().time().timeOutputValue())/Uref_;

        // (PCR:p. 522)
        const vector UBulk
        (
            gSum(UMean()*patch().magSf())
           /(gSum(patch().magSf()) + ROOTVSMALL)
        );

        // Move eddies using bulk velocity
        const scalar deltaT = db().time().deltaTValue();
        convectEddies(UBulk, deltaT);

        // Set mean velocity
        vectorField& U = *this;
        U = UMean;

        // Apply second part of normalisation coefficient
        const scalar c =
            scale_*Foam::pow(10*v0_, m_)/Foam::sqrt(scalar(nEddy_));

        // In parallel, need to collect all eddies that will interact with
        // local faces

        const pointField& Cf = patch().Cf();

        if (singleProc_ || !Pstream::parRun())
        {
            forAll(U, faceI)
            {
                U[faceI] += c*uPrimeEddy(eddies_, Cf[faceI]);
            }
        }
        else
        {
            // Process local eddy contributions
            forAll(U, faceI)
            {
                U[faceI] += c*uPrimeEddy(eddies_, Cf[faceI]);
            }

            // Add contributions from overlapping eddies
            List<List<eddy>> overlappingEddies(Pstream::nProcs());
            calcOverlappingProcEddies(overlappingEddies);

            forAll(overlappingEddies, procI)
            {
                const List<eddy>& eddies = overlappingEddies[procI];

                if (eddies.size())
                {
                    forAll(U, faceI)
                    {
                        U[faceI] += c*uPrimeEddy(eddies, Cf[faceI]);
                    }
                }
            }
        }

        // Re-scale to ensure correct flow rate
        const scalar fCorr =
            gSum((UBulk & patchNormal_)*patch().magSf())
           /gSum(U & -patch().Sf());

        U *= fCorr;

        curTimeIndex_ = db().time().timeIndex();

        if (writeEddies_)
        {
            writeEddyOBJ();
        }

        if (debug)
        {
            Info<< "Magnitude of bulk velocity: " << UBulk << endl;

            label n = eddies_.size();
            Info<< "Number of eddies: " << returnReduce(n, sumOp<label>())
                << endl;

            Info<< "Patch:" << patch().patch().name()
                << " min/max(U):" << gMin(U) << ", " << gMax(U)
                << endl;

            if (db().time().writeTime())
            {
                writeLumleyCoeffs();
            }
        }
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::turbulentDFSEMInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeEntry("delta", delta_);
    os.writeEntryIfDifferent<scalar>("d", 1.0, d_);
    os.writeEntryIfDifferent<scalar>("kappa", 0.41, kappa_);
    os.writeEntryIfDifferent<scalar>("Uref", 1.0, Uref_);
    os.writeEntryIfDifferent<scalar>("Lref", 1.0, Lref_);
    os.writeEntryIfDifferent<scalar>("scale", 1.0, scale_);
    os.writeEntryIfDifferent<scalar>("m", 0.5, m_);
    os.writeEntryIfDifferent<label>("nCellPerEddy", 5, nCellPerEddy_);
    os.writeEntryIfDifferent("writeEddies", false, writeEddies_);
    if (U_)
    {
        U_->writeData(os);
    }
    if (R_)
    {
        R_->writeData(os);
    }
    if (L_)
    {
        L_->writeData(os);
    }
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       turbulentDFSEMInletFvPatchVectorField
   );
}


// ************************************************************************* //
