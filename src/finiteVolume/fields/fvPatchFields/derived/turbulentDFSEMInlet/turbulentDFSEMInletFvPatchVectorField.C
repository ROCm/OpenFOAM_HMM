/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "momentOfInertia.H"
#include "cartesianCS.H"
#include "Fstream.H"
#include "globalIndex.H"

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

        vector offset = patchNormal_*maxSigmaX_;
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

        // Draw lines between points
        // Note: need to order to avoid crossing patch
        //const label nPoint = boundaryPoints.size();
        //
        //forAll(boundaryPoints, i)
        //{
        //    label i1 = i;
        //    label i2 = (i + 1) % nPoint;
        //    os  << "l " << i1 << " " << i2 << nl;
        //}
        //
        //forAll(boundaryPoints, i)
        //{
        //    label i1 = i + nPoint;
        //    label i2 = ((i + 1) % nPoint) + nPoint;
        //    os  << "l " << i1 << " " << i2 << nl;
        //}
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

    // Before interpolation/raw data
    if (interpolateR_)
    {
        fileName valsFile
        (
            fileHandler().filePath
            (
                fileName
                (
                    db().time().path()
                   /db().time().caseConstant()
                   /"boundaryData"
                   /this->patch().name()
                   /"0"
                   /"R"
                )
            )
        );

        autoPtr<ISstream> isPtr
        (
            fileHandler().NewIFstream
            (
                valsFile
            )
        );

        Field<symmTensor> Rexp(isPtr());

        OFstream os(db().time().path()/"lumley_input.out");

        os  << "# xi" << token::TAB << "eta" << endl;

        forAll(Rexp, faceI)
        {
            // Normalised anisotropy tensor
            symmTensor devR = dev(Rexp[faceI]/(tr(Rexp[faceI])));

            // Second tensor invariant
            scalar ii = min(0, invariantII(devR));

            // Third tensor invariant
            scalar iii = invariantIII(devR);

            // xi, eta
            // See Pope - characterization of Reynolds-stress anisotropy
            scalar xi = cbrt(0.5*iii);
            scalar eta = sqrt(-ii/3.0);
            os  << xi << token::TAB << eta << token::TAB
                << ii << token::TAB << iii << endl;
        }
    }

    // After interpolation
    {
        OFstream os(db().time().path()/"lumley_interpolated.out");

        os  << "# xi" << token::TAB << "eta" << endl;

        forAll(R_, faceI)
        {
            // Normalised anisotropy tensor
            symmTensor devR = dev(R_[faceI]/(tr(R_[faceI])));

            // Second tensor invariant
            scalar ii = min(0, invariantII(devR));

            // Third tensor invariant
            scalar iii = invariantIII(devR);

            // xi, eta
            // See Pope - characterization of Reynolds-stress anisotropy
            scalar xi = cbrt(0.5*iii);
            scalar eta = sqrt(-ii/3.0);
            os  << xi << token::TAB << eta << token::TAB
                << ii << token::TAB << iii << endl;
        }
    }
}


const Foam::pointToPointPlanarInterpolation&
Foam::turbulentDFSEMInletFvPatchVectorField::patchMapper() const
{
    // Initialise interpolation (2D planar interpolation by triangulation)
    if (mapperPtr_.empty())
    {
//        vectorGlobalIOField samplePoints
        vectorIOField samplePoints
        (
            IOobject
            (
                "points",
                this->db().time().caseConstant(),
                "boundaryData"/this->patch().name(),
                this->db(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE,
                false
            )
        );

        const fileName samplePointsFile = samplePoints.filePath();

        if (debug)
        {
            InfoInFunction
                << " Read " << samplePoints.size() << " sample points from "
                << samplePointsFile << endl;
        }


        // tbd: run-time selection
        bool nearestOnly =
        (
           !mapMethod_.empty()
         && mapMethod_ != "planarInterpolation"
        );

        // Allocate the interpolator
        mapperPtr_.reset
        (
            new pointToPointPlanarInterpolation
            (
                samplePoints,
                this->patch().patch().faceCentres(),
                perturb_,
                nearestOnly
            )
        );
    }

    return mapperPtr_();
}


void Foam::turbulentDFSEMInletFvPatchVectorField::initialisePatch()
{
    const vectorField nf(patch().nf());

    // Patch normal points into domain
    patchNormal_ = -gAverage(nf);

    // Check that patch is planar
    scalar error = max(magSqr(patchNormal_ + nf));

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

    forAll(sumTriMagSf_, i)
    {
        sumTriMagSf_[i] = 0.0;
    }

    sumTriMagSf_[Pstream::myProcNo() + 1] = sum(triMagSf);

    Pstream::listCombineGather(sumTriMagSf_, maxEqOp<scalar>());
    Pstream::listCombineScatter(sumTriMagSf_);

    for (label i = 1; i < triMagSf.size(); i++)
    {
        triMagSf[i] += triMagSf[i-1];
    }

    // Transfer to persistent storage
    triFace_.transfer(triFace);
    triToFace_.transfer(triToFace);
    triCumulativeMagSf_.transfer(triMagSf);

    // Convert sumTriMagSf_ into cumulative sum of areas per proc
    for (label i = 1; i < sumTriMagSf_.size(); i++)
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

    //const scalarField cellDx(Foam::sqrt(magSf));
    const scalarField cellDx(max(Foam::sqrt(magSf), 2/patch().deltaCoeffs()));

    // Inialise eddy box extents
    forAll(*this, faceI)
    {
        scalar& s = sigmax_[faceI];

        // Length scale in x direction (based on eq. 14)
        s = mag(L_[faceI]);
        s = min(s, kappa_*delta_);

        // Allow eddies to be smaller than the mesh scale as suggested by
        // the reference?
        // s = min(s, nCellPerEddy_*cellDx[faceI]);
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
        scalar areaFraction = rndGen_.globalPosition<scalar>(0, patchArea_);

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
            scalar offset = sumTriMagSf_[procI];
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
        scalar maxAreaLimit = triCumulativeMagSf_.last();
        scalar areaFraction = rndGen_.position<scalar>(0, maxAreaLimit);

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
            label faceI = pos.index();

            // Note: only 1 processor will pick up this face
            if (faceI != -1)
            {
                eddy e
                (
                    faceI,
                    pos.hitPoint(),
                    rndGen_.position<scalar>(-maxSigmaX_, maxSigmaX_),
                    sigmax_[faceI],
                    R_[faceI],
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
            << ": No eddies seeded - please check your set-up" << endl;
    }
}


void Foam::turbulentDFSEMInletFvPatchVectorField::convectEddies
(
    const scalar deltaT
)
{
    // Note: all operations applied to local processor only

    label nRecycled = 0;

    forAll(eddies_, eddyI)
    {
        eddy& e = eddies_[eddyI];
        e.move(deltaT*(UMean_ & patchNormal_));

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
               label faceI = pos.index();

               e = eddy
                    (
                        faceI,
                        pos.hitPoint(),
                        position0 - floor(position0/maxSigmaX_)*maxSigmaX_,
                        sigmax_[faceI],
                        R_[faceI],
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
        Info<< "Patch: " << patch().patch().name() << " recycled "
            << nRecycled << " eddies" << endl;
    }
}


Foam::vector Foam::turbulentDFSEMInletFvPatchVectorField::uDashEddy
(
    const List<eddy>& eddies,
    const point& patchFaceCf
) const
{
    vector uDash(vector::zero);

    forAll(eddies, k)
    {
        const eddy& e = eddies[k];
        uDash += e.uDash(patchFaceCf, patchNormal_);
    }

    return uDash;
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
        point x = e.position(patchNormal_);
        boundBox ebb = e.bounds();
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
            label nRecv = sendSizes[procI][Pstream::myProcNo()];
            constructMap[procI].setSize(nRecv);

            for (label i = 0; i < nRecv; i++)
            {
                constructMap[procI][i] = segmentI++;
            }
        }
    }

    mapDistribute map(segmentI, sendMap.xfer(), constructMap.xfer());

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    for (label domain = 0; domain < Pstream::nProcs(); domain++)
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
    for (label domain = 0; domain < Pstream::nProcs(); domain++)
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
    delta_(0),
    d_(0),
    kappa_(0),

    perturb_(1e-5),
    mapMethod_("planarInterpolation"),
    mapperPtr_(nullptr),
    interpolateR_(false),
    R_(),
    interpolateL_(false),
    L_(),
    interpolateU_(false),
    U_(),
    UMean_(vector::zero),

    patchArea_(-1),
    triFace_(),
    triToFace_(),
    triCumulativeMagSf_(),
    sumTriMagSf_(Pstream::nProcs() + 1, 0.0),

    eddies_(0),
    nCellPerEddy_(5),
    patchNormal_(vector::zero),
    v0_(0),
    rndGen_(Pstream::myProcNo()),
    sigmax_(size(), 0),
    maxSigmaX_(0),
    nEddy_(0),
    curTimeIndex_(-1),
    patchBounds_(boundBox::invertedBox),
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
    delta_(ptf.delta_),
    d_(ptf.d_),
    kappa_(ptf.kappa_),

    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(nullptr),
    interpolateR_(ptf.interpolateR_),
    R_(ptf.R_, mapper),
    interpolateL_(ptf.interpolateL_),
    L_(ptf.L_, mapper),
    interpolateU_(ptf.interpolateU_),
    U_(ptf.U_, mapper),
    UMean_(ptf.UMean_),

    patchArea_(ptf.patchArea_),
    triFace_(ptf.triFace_),
    triToFace_(ptf.triToFace_),
    triCumulativeMagSf_(ptf.triCumulativeMagSf_),
    sumTriMagSf_(ptf.sumTriMagSf_),

    eddies_(ptf.eddies_),
    nCellPerEddy_(ptf.nCellPerEddy_),
    patchNormal_(ptf.patchNormal_),
    v0_(ptf.v0_),
    rndGen_(ptf.rndGen_),
    sigmax_(ptf.sigmax_, mapper),
    maxSigmaX_(ptf.maxSigmaX_),
    nEddy_(0),
    curTimeIndex_(-1),
    patchBounds_(ptf.patchBounds_),
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
    delta_(readScalar(dict.lookup("delta"))),
    d_(dict.lookupOrDefault<scalar>("d", 1)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),

    perturb_(dict.lookupOrDefault<scalar>("perturb", 1e-5)),
    mapMethod_(dict.lookup("mapMethod")),
    mapperPtr_(nullptr),
    interpolateR_(false),
    R_(interpolateOrRead<symmTensor>("R", dict, interpolateR_)),
    interpolateL_(false),
    L_(interpolateOrRead<scalar>("L", dict, interpolateL_)),
    interpolateU_(false),
    U_(interpolateOrRead<vector>("U", dict, interpolateU_)),
    UMean_(vector::zero),

    patchArea_(-1),
    triFace_(),
    triToFace_(),
    triCumulativeMagSf_(),
    sumTriMagSf_(Pstream::nProcs() + 1, 0.0),

    eddies_(),
    nCellPerEddy_(dict.lookupOrDefault<label>("nCellPerEddy", 5)),
    patchNormal_(vector::zero),
    v0_(0),
    rndGen_(0, -1),
    sigmax_(size(), 0),
    maxSigmaX_(0),
    nEddy_(0),
    curTimeIndex_(-1),
    patchBounds_(boundBox::invertedBox),
    singleProc_(false),
    writeEddies_(dict.lookupOrDefault<bool>("writeEddies", false))
{
    eddy::debug = debug;

    // Set UMean as patch area average value
    UMean_ = gSum(U_*patch().magSf())/(gSum(patch().magSf()) + ROOTVSMALL);
}


Foam::turbulentDFSEMInletFvPatchVectorField::
turbulentDFSEMInletFvPatchVectorField
(
    const turbulentDFSEMInletFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    delta_(ptf.delta_),
    d_(ptf.d_),
    kappa_(ptf.kappa_),

    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(nullptr),
    interpolateR_(ptf.interpolateR_),
    R_(ptf.R_),
    interpolateL_(ptf.interpolateL_),
    L_(ptf.L_),
    interpolateU_(ptf.interpolateU_),
    U_(ptf.U_),
    UMean_(ptf.UMean_),

    patchArea_(ptf.patchArea_),
    triFace_(ptf.triFace_),
    triToFace_(ptf.triToFace_),
    triCumulativeMagSf_(ptf.triCumulativeMagSf_),
    sumTriMagSf_(ptf.sumTriMagSf_),

    eddies_(ptf.eddies_),
    nCellPerEddy_(ptf.nCellPerEddy_),
    patchNormal_(ptf.patchNormal_),
    v0_(ptf.v0_),
    rndGen_(ptf.rndGen_),
    sigmax_(ptf.sigmax_),
    maxSigmaX_(ptf.maxSigmaX_),
    nEddy_(0),
    curTimeIndex_(-1),
    patchBounds_(ptf.patchBounds_),
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
    delta_(ptf.delta_),
    d_(ptf.d_),
    kappa_(ptf.kappa_),

    perturb_(ptf.perturb_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(nullptr),
    interpolateR_(ptf.interpolateR_),
    R_(ptf.R_),
    interpolateL_(ptf.interpolateL_),
    L_(ptf.L_),
    interpolateU_(ptf.interpolateU_),
    U_(ptf.U_),
    UMean_(ptf.UMean_),

    patchArea_(ptf.patchArea_),
    triFace_(ptf.triFace_),
    triToFace_(ptf.triToFace_),
    triCumulativeMagSf_(ptf.triCumulativeMagSf_),
    sumTriMagSf_(ptf.sumTriMagSf_),

    eddies_(ptf.eddies_),
    nCellPerEddy_(ptf.nCellPerEddy_),
    patchNormal_(ptf.patchNormal_),
    v0_(ptf.v0_),
    rndGen_(ptf.rndGen_),
    sigmax_(ptf.sigmax_),
    maxSigmaX_(ptf.maxSigmaX_),
    nEddy_(0),
    curTimeIndex_(-1),
    patchBounds_(ptf.patchBounds_),
    singleProc_(ptf.singleProc_),
    writeEddies_(ptf.writeEddies_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turbulentDFSEMInletFvPatchVectorField::~
turbulentDFSEMInletFvPatchVectorField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::turbulentDFSEMInletFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<vector>::autoMap(m);

    // Clear interpolator
    mapperPtr_.clear();
    R_.autoMap(m);
    L_.autoMap(m);
    U_.autoMap(m);

    sigmax_.autoMap(m);
}


void Foam::turbulentDFSEMInletFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<vector>::rmap(ptf, addr);

    const turbulentDFSEMInletFvPatchVectorField& dfsemptf =
        refCast<const turbulentDFSEMInletFvPatchVectorField>(ptf);

    R_.rmap(dfsemptf.R_, addr);
    L_.rmap(dfsemptf.L_, addr);
    U_.rmap(dfsemptf.U_, addr);

    // Clear interpolator
    mapperPtr_.clear();

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
        if (debug)
        {
            label n = eddies_.size();
            Info<< "Number of eddies: " << returnReduce(n, sumOp<label>())
                << endl;
        }

        const scalar deltaT = db().time().deltaTValue();

        // Move eddies using mean velocity
        convectEddies(deltaT);

        // Set velocity
        vectorField& U = *this;
        //U = UMean_;
        U = U_;

        const pointField& Cf = patch().Cf();

        // Apply second part of normalisation coefficient
        // Note: factor of 2 required to match reference stresses?
        const scalar FACTOR = 2;
        const scalar c = FACTOR*Foam::sqrt(10*v0_)/Foam::sqrt(scalar(nEddy_));

        // In parallel, need to collect all eddies that will interact with
        // local faces

        if (singleProc_ || !Pstream::parRun())
        {
            forAll(U, faceI)
            {
                U[faceI] += c*uDashEddy(eddies_, Cf[faceI]);
            }
        }
        else
        {
            // Process local eddy contributions
            forAll(U, faceI)
            {
                U[faceI] += c*uDashEddy(eddies_, Cf[faceI]);
            }

            // Add contributions from overlapping eddies
            List<List<eddy>> overlappingEddies(Pstream::nProcs());
            calcOverlappingProcEddies(overlappingEddies);

            forAll(overlappingEddies, procI)
            {
                const List<eddy>& eddies = overlappingEddies[procI];

                if (eddies.size())
                {
                    //Pout<< "Applying " << eddies.size()
                    //    << " eddies from processor " << procI << endl;

                    forAll(U, faceI)
                    {
                        U[faceI] += c*uDashEddy(eddies, Cf[faceI]);
                    }
                }
            }
        }

        // Re-scale to ensure correct flow rate
        scalar fCorr =
            gSum((UMean_ & patchNormal_)*patch().magSf())
           /gSum(U & -patch().Sf());

        U *= fCorr;

        if (debug)
        {
            Info<< "Patch:" << patch().patch().name()
                << " min/max(U):" << gMin(U) << ", " << gMax(U) << endl;
        }

        curTimeIndex_ = db().time().timeIndex();

        if (writeEddies_)
        {
            writeEddyOBJ();
        }

        if (debug && db().time().outputTime())
        {
            writeLumleyCoeffs();
        }
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::turbulentDFSEMInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    writeEntry("value", os);
    os.writeEntry("delta", delta_);
    os.writeEntryIfDifferent<scalar>("d", 1.0, d_);
    os.writeEntryIfDifferent<scalar>("kappa", 0.41, kappa_);
    os.writeEntryIfDifferent<scalar>("perturb", 1e-5, perturb_);
    os.writeEntryIfDifferent<label>("nCellPerEddy", 5, nCellPerEddy_);
    os.writeEntryIfDifferent("writeEddies", false, writeEddies_);

    if (!interpolateR_)
    {
        R_.writeEntry("R", os);
    }

    if (!interpolateL_)
    {
        L_.writeEntry("L", os);
    }

    if (!interpolateU_)
    {
        U_.writeEntry("U", os);
    }

    if (!mapMethod_.empty())
    {
        os.writeEntryIfDifferent<word>
        (
            "mapMethod",
            "planarInterpolation",
            mapMethod_
        );
    }
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
