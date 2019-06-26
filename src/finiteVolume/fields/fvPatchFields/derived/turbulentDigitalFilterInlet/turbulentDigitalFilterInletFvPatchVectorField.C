/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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

#include "turbulentDigitalFilterInletFvPatchVectorField.H"
#include "volFields.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::turbulentDigitalFilterInletFvPatchVectorField::variantType
>
Foam::turbulentDigitalFilterInletFvPatchVectorField::variantNames
({
    { variantType::DIGITAL_FILTER, "digitalFilter" },
    { variantType::DIGITAL_FILTER, "DFM" },
    { variantType::FORWARD_STEPWISE, "forwardStepwise" },
    { variantType::FORWARD_STEPWISE, "reducedDigitalFilter" },
    { variantType::FORWARD_STEPWISE, "FSM" }
});

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::pointToPointPlanarInterpolation&
Foam::turbulentDigitalFilterInletFvPatchVectorField::patchMapper() const
{
    // Initialise interpolation (2D planar interpolation by triangulation)
    if (mapperPtr_.empty())
    {
        // Reread values and interpolate
        fileName samplePointsFile
        (
            this->db().time().path()
           /this->db().time().caseConstant()
           /"boundaryData"
           /this->patch().name()
           /"points"
        );

        pointField samplePoints((IFstream(samplePointsFile)()));

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

    return *mapperPtr_;
}


Foam::List<Foam::Pair<Foam::label>>
Foam::turbulentDigitalFilterInletFvPatchVectorField::patchIndexPairs()
{
    // Get patch normal direction into the domain
    const vector nf(computePatchNormal());

    // Find the second local coordinate direction
    direction minCmpt = -1;
    scalar minMag = VGREAT;
    for (direction cmpt = 0; cmpt < pTraits<vector>::nComponents; ++cmpt)
    {
        scalar s = mag(nf[cmpt]);
        if (s < minMag)
        {
            minMag = s;
            minCmpt = cmpt;
        }
    }

    // Create the second local coordinate direction
    vector e2(Zero);
    e2[minCmpt] = 1.0;

    // Remove normal component
    e2 -= (nf&e2)*nf;

    // Create the local coordinate system
    coordSystem::cartesian cs
    (
        Zero,   // origin
        nf,     // normal
        e2      // 0-axis
    );

    // Convert patch points into local coordinate system
    const pointField localPos
    (
        cs.localPosition
        (
            pointField
            (
                patch().patch().points(),
                patch().patch().meshPoints()
            )
        )
    );

    const boundBox bb(localPos);

    // Normalise to (i, j) coordinates
    const Vector2D<label> n(planeDivisions_.first(), planeDivisions_.second());
    invDelta_[0] = n[0]/bb.span()[0];
    invDelta_[1] = n[1]/bb.span()[1];
    const point localMinPt(bb.min());

    // Compute virtual-actual patch index pairs
    List<Pair<label>> indexPairs(this->size(), Pair<label>(Zero, Zero));

    // Virtual turbulence plane indices
    label j = 0;
    label k = 0;

    forAll(*this, facei)
    {
        const scalar& centre0 = localPos[facei][0];
        const scalar& centre1 = localPos[facei][1];

        j = label((centre0 - localMinPt[0])*invDelta_[0]);
        k = label((centre1 - localMinPt[1])*invDelta_[1]);

        indexPairs[facei] = Pair<label>(facei, k*n[0] + j);
    }

    return indexPairs;
}


void Foam::turbulentDigitalFilterInletFvPatchVectorField::
checkRTensorRealisable() const
{
    const vectorField& faceCentres = this->patch().patch().faceCentres();

    forAll(R_, facei)
    {
        if (R_[facei].xx() <= 0)
        {
            FatalErrorInFunction
                << "Reynolds stress tensor component Rxx cannot be negative"
                   "or zero, where Rxx = " << R_[facei].xx() << " at the face "
                   "centre = " << faceCentres[facei] << exit(FatalError);
        }

        if (R_[facei].yy() < 0 || R_[facei].zz() < 0)
        {
            FatalErrorInFunction
                << "Reynolds stress tensor components Ryy or Rzz cannot be"
                << "negative where Ryy = " << R_[facei].yy() << ", and Rzz = "
                << R_[facei].zz() << " at the face centre = "
                << faceCentres[facei] << exit(FatalError);
        }

        scalar term0 = R_[facei].xx()*R_[facei].yy() - sqr(R_[facei].xy());

        if (term0 <= 0)
        {
            FatalErrorInFunction
                << "Reynolds stress tensor component group, Rxx*Ryy - Rxy^2"
                << "cannot be negative or zero at the face centre = "
                << faceCentres[facei] << exit(FatalError);
        }

        scalar term1 = R_[facei].zz() - sqr(R_[facei].xz())/R_[facei].xx();
        scalar term2 =
            sqr(R_[facei].yz() - R_[facei].xy()*R_[facei].xz()
           /(R_[facei].xx()*term0));
        scalar term3 = term1 - term2;

        if (term3 < 0)
        {
            FatalErrorInFunction
                << "Reynolds stress tensor component group,"
                << "Rzz - Rxz^2/Rxx - (Ryz - Rxy*Rxz/(Rxx*(Rxx*Ryy - Rxy^2)))^2"
                << "cannot be negative at the face centre = "
                << faceCentres[facei] << exit(FatalError);
        }
    }

    #ifdef FULLDEBUG
    Info<< "Ends: checkRTensorRealisable()."
        << " Reynolds tensor (on patch) is consistent." << nl;
    #endif
}


Foam::symmTensorField Foam::turbulentDigitalFilterInletFvPatchVectorField::
computeLundWuSquires() const
{
    checkRTensorRealisable();

    symmTensorField LundWuSquires(symmTensorField(R_.size()));

    forAll(R_, facei)
    {
        const symmTensor& R = R_[facei];
        symmTensor& lws = LundWuSquires[facei];

        // (Klein et al., 2003, Eq. 5)
        lws.xx() = Foam::sqrt(R.xx());
        lws.xy() = R.xy()/lws.xx();
        lws.xz() = R.xz()/lws.xx();
        lws.yy() = Foam::sqrt(R.yy() - sqr(lws.xy()));
        lws.yz() = (R.yz() - lws.xy()*lws.xz())/lws.yy();
        lws.zz() = Foam::sqrt(R.zz() - sqr(lws.xz()) - sqr(lws.yz()));
    }

    #ifdef FULLDEBUG
    Info<< "Ends: computeLundWuSquires()." << nl;
    #endif

    return LundWuSquires;
}


Foam::vector Foam::turbulentDigitalFilterInletFvPatchVectorField::
computePatchNormal() const
{
    vector patchNormal(-gAverage(patch().nf()));
    return patchNormal.normalise();
}


Foam::scalar Foam::turbulentDigitalFilterInletFvPatchVectorField::
computeInitialFlowRate() const
{
    const vector patchNormal(computePatchNormal());
    return gSum((UMean_ & patchNormal)*patch().magSf());
}


void Foam::turbulentDigitalFilterInletFvPatchVectorField::convertToTimeScale
(
    tensor& L
) const
{
    if (isTaylorHypot_)
    {
        forAll(L.x(), i)
        {
            L[i] /= patchNormalSpeed_;
        }

    #ifdef FULLDEBUG
    Info<< "Ends: convertToTimeScale()."
        << "Streamwise integral length scales are converted to time scales via"
        << "Taylor's frozen turbulence hypothesis" << nl;
    #endif
    }
}


Foam::tensor Foam::turbulentDigitalFilterInletFvPatchVectorField::
convertScalesToGridUnits
(
    const tensor& L
) const
{
    tensor Ls(L);
    convertToTimeScale(Ls);
    const scalar invDeltaT = 1.0/db().time().deltaTValue();

    Ls.row(0, Ls.x()*invDeltaT);
    Ls.row(1, Ls.y()*invDelta_[0]);
    Ls.row(2, Ls.z()*invDelta_[1]);

    #ifdef FULLDEBUG
    Info<< "Ends: convertScalesToGridUnits()."
        << "Units of input length scales are converted from metres to"
        << "virtual-patch cell size/time-step" << nl;
    #endif

    return Ls;
}


Foam::List<Foam::label> Foam::turbulentDigitalFilterInletFvPatchVectorField::
initLenRandomBox() const
{
    label initValue = 0;
    label rangeModifier = 0;

    if (variant_ == variantType::FORWARD_STEPWISE)
    {
        // Initialise with 1 since x-dir possess 1 node with this variant
        initValue = pTraits<label>::nComponents;
        rangeModifier = pTraits<vector>::nComponents;
    }

    List<label> lenRandomBox(pTraits<tensor>::nComponents, initValue);
    Vector<label> lenGrid
    (
        pTraits<label>::nComponents,
        planeDivisions_.first(),
        planeDivisions_.second()
    );

    // Initialise: Main convenience factor, lenRandomBox_
    for
    (
        const label& i
      : labelRange(rangeModifier, pTraits<tensor>::nComponents - rangeModifier)
    )
    {
        // Slicing index
        const label sliceI = label(i/pTraits<vector>::nComponents);

        // Refer to 'computeFilterCoeffs()'
        const label n = ceil(L_[i]);
        const label twiceN = 4*n;

        // Initialise: Random-number set sizes
        lenRandomBox[i] = lenGrid[sliceI] + twiceN;
    }

    return lenRandomBox;
}


Foam::List<Foam::label> Foam::turbulentDigitalFilterInletFvPatchVectorField::
initBoxFactors2D() const
{
    List<label> boxFactors2D(pTraits<vector>::nComponents);

    for (direction dir = 0; dir < pTraits<vector>::nComponents; ++dir)
    {
        boxFactors2D[dir] =
            lenRandomBox_[pTraits<vector>::nComponents + dir]
           *lenRandomBox_[pTraits<symmTensor>::nComponents + dir];
    }

    return boxFactors2D;
}


Foam::List<Foam::label> Foam::turbulentDigitalFilterInletFvPatchVectorField::
initBoxFactors3D() const
{
    List<label> boxFactors3D(pTraits<vector>::nComponents);

    for (direction dir = 0; dir < pTraits<vector>::nComponents; ++dir)
    {
        boxFactors3D[dir] = randomBoxFactors2D_[dir]*lenRandomBox_[dir];
    }

    return boxFactors3D;
}


Foam::List<Foam::label> Foam::turbulentDigitalFilterInletFvPatchVectorField::
initBoxPlaneFactors() const
{
    List<label> planeFactors(pTraits<vector>::nComponents);

    for (direction dir = 0; dir < pTraits<vector>::nComponents; ++dir)
    {
        planeFactors[dir] =
            randomBoxFactors2D_[dir]*(lenRandomBox_[dir] - pTraits<label>::one);
    }

    return planeFactors;
}


Foam::List<Foam::List<Foam::scalar>>
Foam::turbulentDigitalFilterInletFvPatchVectorField::fillRandomBox()
{
    List<List<scalar>> randomBox(pTraits<vector>::nComponents, List<scalar>());

    // Initialise: Remaining convenience factors for (e1 e2 e3)
    for (direction dir = 0; dir < pTraits<vector>::nComponents; ++dir)
    {
        // Initialise: Random-box content with random-number sets
        randomBox[dir] =
            generateRandomSet<List<scalar>, scalar>(randomBoxFactors3D_[dir]);
    }

    #ifdef FULLDEBUG
    Info<< "Ends: fillRandomBox()."
        << "Random-number box is created." << nl;
    #endif

    return randomBox;
}


Foam::List<Foam::List<Foam::scalar>>
Foam::turbulentDigitalFilterInletFvPatchVectorField::computeFilterCoeffs() const
{
    List<List<scalar>> filterCoeffs
    (
        pTraits<tensor>::nComponents,
        List<scalar>(1, 1.0)
    );

    label initValue = 0;

    if(variant_ == variantType::FORWARD_STEPWISE)
    {
        initValue = pTraits<vector>::nComponents;
    }

    for (direction dir = initValue; dir < pTraits<tensor>::nComponents; ++dir)
    {
        // The smallest filter width larger than length scale
        // (Klein et al., 2003, 'n' in Eq. 15)
        const label n  = ceil(L_[dir]);

        // Full filter-width (except mid-zero) according to the condition
        // (Klein et al., 2003, Eq. 15 whereat N is minimum =2n)
        const label twiceN = 4*n;

        // Initialise filter-coeffs containers with full filter-width size
        // Extra elem is mid-zero within [-N, N]
        filterCoeffs[dir] = List<scalar>(twiceN + 1, Zero);

        // First element: -N within [-N, N]
        const scalar initElem = -2*scalar(n);

        // Container initialised with [-N, N]
        // (Klein et al., 2003, p. 658, Item-b)
        std::iota
        (
            filterCoeffs[dir].begin(),
            filterCoeffs[dir].end(),
            initElem
        );

        // Compute filter-coeffs
        // (Klein et al., 2003, Eq. 14 (Gaussian))
        // (Bercin et al., 2018, Fig. 9 (Exponential))
        List<scalar> fTemp(filterCoeffs[dir]);
        scalar fSum = 0;
        const scalar nSqr = n*n;

        if (isGaussian_)
        {
            fTemp = sqr(exp(modelConst_*sqr(fTemp)/nSqr));
            fSum = Foam::sqrt(sum(fTemp));

            filterCoeffs[dir] =
                exp(modelConst_*sqr(filterCoeffs[dir])/nSqr)/fSum;
        }
        else
        {
            fTemp = sqr(exp(modelConst_*mag(fTemp)/n));
            fSum = Foam::sqrt(sum(fTemp));

            filterCoeffs[dir] = exp(modelConst_*mag(filterCoeffs[dir])/n)/fSum;
        }
    }

    #ifdef FULLDEBUG
    Info<< "Ends: computeFilterCoeffs()."
        << " Filter coefficients are computed." << nl;
    #endif

    return filterCoeffs;
}


void Foam::turbulentDigitalFilterInletFvPatchVectorField::rndShiftRefill()
{
    forAll(randomBox_, dir)
    {
        List<scalar>& r = randomBox_[dir];

        // Shift forward from the back to the front / First Out
        inplaceRotateList(r, randomBoxFactors2D_[dir]);

        // Refill the back with a new random-number set / First In
        for (label i = 0; i < randomBoxFactors2D_[dir]; ++i)
        {
            r[i] = rndGen_.GaussNormal<scalar>();
        }
    }
}


void Foam::turbulentDigitalFilterInletFvPatchVectorField::mapFilteredRandomBox
(
    vectorField& U
)
{
    for (const auto& x : indexPairs_)
    {
        const label xf = x.first();
        const label xs = x.second();
        U[xf][0] = filteredRandomBox_[0][xs];
        U[xf][1] = filteredRandomBox_[1][xs];
        U[xf][2] = filteredRandomBox_[2][xs];
    }
}


void Foam::turbulentDigitalFilterInletFvPatchVectorField::embedOnePointCorrs
(
    vectorField& U
) const
{
    forAll(LundWuSquires_, facei)
    {
        vector& Us = U[facei];
        const symmTensor& lws = LundWuSquires_[facei];

        // (Klein et al. p. 658, Item-e)
        Us.z() = Us.x()*lws.xz() + Us.y()*lws.yz() + Us.z()*lws.zz();
        Us.y() = Us.x()*lws.xy() + Us.y()*lws.yy();
        Us.x() = Us.x()*lws.xx();
    }
}


void Foam::turbulentDigitalFilterInletFvPatchVectorField::embedMeanVelocity
(
    vectorField& U
) const
{
    U += UMean_;
}


void Foam::turbulentDigitalFilterInletFvPatchVectorField::correctFlowRate
(
    vectorField& U
) const
{
    U *= (initialFlowRate_/gSum(U & -patch().Sf()));
}


void Foam::turbulentDigitalFilterInletFvPatchVectorField::embedTwoPointCorrs()
{
    for (direction dir = 0; dir < pTraits<vector>::nComponents; ++dir)
    {
        List<scalar>& in = randomBox_[dir];
        List<scalar>& out = filteredRandomBox_[dir];
        const List<scalar>& filter1 = filterCoeffs_[dir];
        const List<scalar>& filter2 = filterCoeffs_[3 + dir];
        const List<scalar>& filter3 = filterCoeffs_[6 + dir];

        const label sz1 = lenRandomBox_[dir];
        const label sz2 = lenRandomBox_[3 + dir];
        const label sz3 = lenRandomBox_[6 + dir];
        const label szfilter1 = filterCoeffs_[dir].size();
        const label szfilter2 = filterCoeffs_[3 + dir].size();
        const label szfilter3 = filterCoeffs_[6 + dir].size();
        const label sz23 = randomBoxFactors2D_[dir];
        const label sz123 = randomBoxFactors3D_[dir];
        const label validSlice2 = sz2 - (szfilter2 - label(1));
        const label validSlice3 = sz3 - (szfilter3 - label(1));

        // Convolution summation - Along 1st direction
        scalarField tmp(sz123);
        label filterCentre = label(szfilter2/label(2));
        label endIndex = sz2 - filterCentre;
        label i0 = 0;
        label i1 = 0;
        label i2 = 0;

        for (label i = 0; i < sz1; ++i)
        {
            for (label j = 0; j < sz3; ++j)
            {
                i1 += filterCentre;

                for (label k = filterCentre; k < endIndex; ++k)
                {
                    tmp[i1] = 0.0;
                    label q = 0;

                    for (label p = szfilter2 - 1; p >= 0; --p, ++q)
                    {
                        tmp[i1] += in[i0 + q]*filter2[p];
                    }
                    ++i0;
                    ++i1;

                }
                i0 += (filterCentre + filterCentre);
                i1 += filterCentre;
            }
        }

        // Convolution summation - Along 2nd direction
        scalarField tmp2(tmp);
        filterCentre = label(szfilter3/label(2));
        endIndex = sz3 - filterCentre;
        i1 = 0;
        i2 = 0;

        for (label i = 0; i < sz1; ++i)
        {
           label sl = i*sz23;

            for (label j = 0; j < sz2; ++j)
            {
                i1 = j + sl;
                i2 = i1;

                for (label k = 0; k < endIndex - filterCentre; ++k)
                {
                    tmp[i1] = 0.0;
                    label q = 0;

                    for (label p = szfilter3 - 1; p >= 0; --p, ++q)
                    {
                        tmp[i1] += tmp2[i2 + q*sz2]*filter3[p];
                    }
                    i1 += sz2;
                    i2 += sz2;

                }
                i1 += (sz2 + filterCentre);
                i2 += (sz2 + filterCentre);
            }
        }

        // Convolution summation - Along 3rd direction
        filterCentre = label(szfilter1/label(2));
        endIndex = sz1 - filterCentre;
        i1 = (szfilter2 - label(1))/label(2);
        i2 = (szfilter2 - label(1))/label(2);
        label i3 = 0;

        for (label i = 0; i < validSlice3; ++i)
        {
            for (label j = 0; j < validSlice2; ++j)
            {
                scalar sum = 0.0;
                i1 = i2 + j;

                for (label k = szfilter1 - 1; k >= 0; --k)
                {
                    sum += tmp[i1]*filter1[k];
                    i1 += sz23;
                }
                out[i3] = sum;
                ++i3;

            }
            i2 += sz2;
        }
    }
}


void Foam::turbulentDigitalFilterInletFvPatchVectorField::computeDFM
(
    vectorField& U
)
{
    #ifdef FULLDEBUG
    Info<< "Starts: computeDFM()" << nl;
    #endif

    if (Pstream::master())
    {
        embedTwoPointCorrs();
        rndShiftRefill();
    }

    Pstream::scatter(filteredRandomBox_);

    mapFilteredRandomBox(U);

    embedOnePointCorrs(U);

    embedMeanVelocity(U);

    if (isCorrectedFlowRate_)
    {
        correctFlowRate(U);
    }

    #ifdef FULLDEBUG
    Info<< "Ends: computeDFM()" << nl;
    #endif
}


void Foam::turbulentDigitalFilterInletFvPatchVectorField::computeReducedDFM
(
    vectorField& U
)
{
    #ifdef FULLDEBUG
    Info<< "Starts: computeReducedDFM()" << nl;
    #endif

    if (Pstream::master())
    {
        embedTwoPointCorrs();
        rndShiftRefill();
    }

    Pstream::scatter(filteredRandomBox_);

    mapFilteredRandomBox(U);

    computeFSM(U);

    embedOnePointCorrs(U);

    embedMeanVelocity(U);

    if (isCorrectedFlowRate_)
    {
        correctFlowRate(U);
    }

    #ifdef FULLDEBUG
    Info<< "Ends: computeReducedDFM()" << nl;
    #endif
}


Foam::List<Foam::scalar> Foam::turbulentDigitalFilterInletFvPatchVectorField::
computeConstList1FSM() const
{
    List<scalar> constList1FSM(pTraits<vector>::nComponents);

    forAll(L_.x(), i)
    {
        constList1FSM[i] = exp(const1FSM_/L_[i]);
    }

    return constList1FSM;
}


Foam::List<Foam::scalar> Foam::turbulentDigitalFilterInletFvPatchVectorField::
computeConstList2FSM() const
{
    List<scalar> constList2FSM(pTraits<vector>::nComponents);

    forAll(L_.x(), i)
    {
        constList2FSM[i] = Foam::sqrt(1.0 - exp(const2FSM_/L_[i]));
    }

    return constList2FSM;
}


void Foam::turbulentDigitalFilterInletFvPatchVectorField::computeFSM
(
    vectorField& U
)
{
    for (direction dir = 0; dir < pTraits<vector>::nComponents; ++dir)
    {
        U.component(dir) = U0_.component(dir)*constList1FSM_[dir]
                          +U.component(dir)*constList2FSM_[dir];
    }

    U0_ = U;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentDigitalFilterInletFvPatchVectorField::
turbulentDigitalFilterInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    mapperPtr_(nullptr),
    variant_(variantType::DIGITAL_FILTER),
    isGaussian_(true),
    isFixedSeed_(true),
    isContinuous_(false),
    isCorrectedFlowRate_(true),
    interpolateR_(false),
    interpolateUMean_(false),
    isInsideMesh_(false),
    isTaylorHypot_(true),
    mapMethod_("nearestCell"),
    curTimeIndex_(-1),
    tiny_(1e-8),
    patchNormalSpeed_(Zero),
    modelConst_(-0.5*constant::mathematical::pi),
    perturb_(1e-5),
    initialFlowRate_(pTraits<scalar>::one),
    rndGen_(1234),
    planeDivisions_(Zero, Zero),
    invDelta_(Zero),
    indexPairs_(Zero),
    R_(Zero),
    LundWuSquires_(Zero),
    UMean_(Zero),
    Lbak_(Zero),
    L_(Zero),
    const1FSM_(Zero),
    const2FSM_(Zero),
    constList1FSM_(Zero),
    constList2FSM_(Zero),
    lenRandomBox_(Zero),
    randomBoxFactors2D_(Zero),
    randomBoxFactors3D_(Zero),
    iNextToLastPlane_(Zero),
    randomBox_(Zero),
    filterCoeffs_(Zero),
    filteredRandomBox_(Zero),
    U0_(Zero),
    computeVariant(nullptr)
{}


Foam::turbulentDigitalFilterInletFvPatchVectorField::
turbulentDigitalFilterInletFvPatchVectorField
(
    const turbulentDigitalFilterInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    mapperPtr_(nullptr),
    variant_(ptf.variant_),
    isGaussian_(ptf.isGaussian_),
    isFixedSeed_(ptf.isFixedSeed_),
    isContinuous_(ptf.isContinuous_),
    isCorrectedFlowRate_(ptf.isCorrectedFlowRate_),
    interpolateR_(ptf.interpolateR_),
    interpolateUMean_(ptf.interpolateUMean_),
    isInsideMesh_(ptf.isInsideMesh_),
    isTaylorHypot_(ptf.isTaylorHypot_),
    mapMethod_(ptf.mapMethod_),
    curTimeIndex_(-1),
    tiny_(ptf.tiny_),
    patchNormalSpeed_(ptf.patchNormalSpeed_),
    modelConst_(ptf.modelConst_),
    perturb_(ptf.perturb_),
    initialFlowRate_(ptf.initialFlowRate_),
    rndGen_(ptf.rndGen_),
    planeDivisions_(ptf.planeDivisions_),
    invDelta_(ptf.invDelta_),
    indexPairs_(ptf.indexPairs_),
    R_(ptf.R_),
    LundWuSquires_(ptf.LundWuSquires_),
    UMean_(ptf.UMean_),
    Lbak_(ptf.Lbak_),
    L_(ptf.L_),
    const1FSM_(ptf.const1FSM_),
    const2FSM_(ptf.const2FSM_),
    constList1FSM_(ptf.constList1FSM_),
    constList2FSM_(ptf.constList2FSM_),
    lenRandomBox_(ptf.lenRandomBox_),
    randomBoxFactors2D_(ptf.randomBoxFactors2D_),
    randomBoxFactors3D_(ptf.randomBoxFactors3D_),
    iNextToLastPlane_(ptf.iNextToLastPlane_),
    randomBox_(ptf.randomBox_),
    filterCoeffs_(ptf.filterCoeffs_),
    filteredRandomBox_(ptf.filteredRandomBox_),
    U0_(ptf.U0_),
    computeVariant(ptf.computeVariant)
{}


Foam::turbulentDigitalFilterInletFvPatchVectorField::
turbulentDigitalFilterInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    mapperPtr_(nullptr),
    variant_
    (
        variantNames.getOrDefault
        (
            "variant",
            dict,
            variantType::DIGITAL_FILTER
        )
    ),
    isGaussian_
    (
        (variant_ == variantType::FORWARD_STEPWISE)
      ? false
      : dict.getOrDefault("isGaussian", true)
    ),
    isFixedSeed_(dict.getOrDefault("isFixedSeed", true)),
    isContinuous_(dict.getOrDefault("isContinuous", false)),
    isCorrectedFlowRate_(dict.getOrDefault("isCorrectedFlowRate", true)),
    interpolateR_(dict.getOrDefault("interpolateR", false)),
    interpolateUMean_(dict.getOrDefault("interpolateUMean", false)),
    isInsideMesh_(dict.getOrDefault("isInsideMesh", false)),
    isTaylorHypot_(dict.getOrDefault("isTaylorHypot", true)),
    mapMethod_(dict.getOrDefault<word>("mapMethod", "nearestCell")),
    curTimeIndex_(-1),
    tiny_(dict.getOrDefault<scalar>("threshold", 1e-8)),
    patchNormalSpeed_
    (
        dict.getCheck<scalar>
        (
            "patchNormalSpeed",
            scalarMinMax::ge(tiny_)
        )
    ),
    modelConst_
    (
        dict.getOrDefault<scalar>
        (
            "modelConst",
           -0.5*constant::mathematical::pi
        )
    ),
    perturb_(dict.getOrDefault<scalar>("perturb", 1e-5)),
    initialFlowRate_(pTraits<scalar>::one),
    rndGen_
    (
        isFixedSeed_
      ? 1234
      : time(0)
    ),
    planeDivisions_
    (
        dict.getCheck<Tuple2<label, label>>
        (
            "planeDivisions",
            [&](const Tuple2<label, label>& len)
            {
                return tiny_ < min(len.first(), len.second()) ? true : false;
            }
        )
    ),
    invDelta_(),
    indexPairs_(patchIndexPairs()),
    R_(interpolateOrRead<symmTensor>("R", dict, interpolateR_)),
    LundWuSquires_(computeLundWuSquires()),
    UMean_(interpolateOrRead<vector>("UMean", dict, interpolateUMean_)),
    Lbak_
    (
        dict.getCheck<tensor>
        (
            "L",
            [&](const tensor& l){return tiny_ < cmptMin(l) ? true : false;}
        )
    ),
    L_(convertScalesToGridUnits(Lbak_)),
    const1FSM_
    (
        (variant_ == variantType::FORWARD_STEPWISE)
      ? dict.getOrDefault<scalar>
        (
            "const1FSM",
           -0.25*constant::mathematical::pi
        )
      : Zero
    ),
    const2FSM_
    (
        (variant_ == variantType::FORWARD_STEPWISE)
      ? dict.getOrDefault<scalar>
        (
            "const2FSM",
           -0.5*constant::mathematical::pi
        )
      : Zero
    ),
    constList1FSM_
    (
        (variant_ == variantType::FORWARD_STEPWISE)
      ? computeConstList1FSM()
      : List<scalar>()
    ),
    constList2FSM_
    (
        (variant_ == variantType::FORWARD_STEPWISE)
      ? computeConstList2FSM()
      : List<scalar>()
    ),
    lenRandomBox_(initLenRandomBox()),
    randomBoxFactors2D_(initBoxFactors2D()),
    randomBoxFactors3D_(initBoxFactors3D()),
    iNextToLastPlane_(initBoxPlaneFactors()),
    randomBox_
    (
        (isContinuous_ && Pstream::master())
      ? dict.getOrDefault<List<List<scalar>>>
        (
            "randomBox",
            fillRandomBox() // First time-step
        )
      :
        fillRandomBox()
    ),
    filterCoeffs_(computeFilterCoeffs()),
    filteredRandomBox_
    (
        pTraits<vector>::nComponents,
        List<scalar>(planeDivisions_.first()*planeDivisions_.second(), Zero)
    ),
    U0_
    (
        (variant_ == variantType::FORWARD_STEPWISE)
      ? generateRandomSet<vectorField, vector>(patch().size())
      : vectorField()
    ),
    computeVariant
    (
        (variant_ == variantType::FORWARD_STEPWISE)
      ? &turbulentDigitalFilterInletFvPatchVectorField::computeReducedDFM
      : &turbulentDigitalFilterInletFvPatchVectorField::computeDFM
    )
{
    if (isCorrectedFlowRate_)
    {
        initialFlowRate_ = computeInitialFlowRate();
    }

    // Check if varying or fixed time-step computation
    if (db().time().isAdjustTimeStep())
    {
        WarningInFunction
            << "Varying time-step computations are not fully supported"
            << " for the moment."<< nl << nl;
    }

    #ifdef FULLDEBUG
    Info<< "Ends: Resource acquisition/initialisation for the"
        << " synthetic turbulence boundary condition." << nl;
    #endif
}


Foam::turbulentDigitalFilterInletFvPatchVectorField::
turbulentDigitalFilterInletFvPatchVectorField
(
    const turbulentDigitalFilterInletFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    mapperPtr_(nullptr),
    variant_(ptf.variant_),
    isGaussian_(ptf.isGaussian_),
    isFixedSeed_(ptf.isFixedSeed_),
    isContinuous_(ptf.isContinuous_),
    isCorrectedFlowRate_(ptf.isCorrectedFlowRate_),
    interpolateR_(ptf.interpolateR_),
    interpolateUMean_(ptf.interpolateUMean_),
    isInsideMesh_(ptf.isInsideMesh_),
    isTaylorHypot_(ptf.isTaylorHypot_),
    mapMethod_(ptf.mapMethod_),
    curTimeIndex_(ptf.curTimeIndex_),
    tiny_(ptf.tiny_),
    patchNormalSpeed_(ptf.patchNormalSpeed_),
    modelConst_(ptf.modelConst_),
    perturb_(ptf.perturb_),
    initialFlowRate_(ptf.initialFlowRate_),
    rndGen_(ptf.rndGen_),
    planeDivisions_(ptf.planeDivisions_),
    invDelta_(ptf.invDelta_),
    indexPairs_(ptf.indexPairs_),
    R_(ptf.R_),
    LundWuSquires_(ptf.LundWuSquires_),
    UMean_(ptf.UMean_),
    Lbak_(ptf.Lbak_),
    L_(ptf.L_),
    const1FSM_(ptf.const1FSM_),
    const2FSM_(ptf.const2FSM_),
    constList1FSM_(ptf.constList1FSM_),
    constList2FSM_(ptf.constList2FSM_),
    lenRandomBox_(ptf.lenRandomBox_),
    randomBoxFactors2D_(ptf.randomBoxFactors2D_),
    randomBoxFactors3D_(ptf.randomBoxFactors3D_),
    iNextToLastPlane_(ptf.iNextToLastPlane_),
    randomBox_(ptf.randomBox_),
    filterCoeffs_(ptf.filterCoeffs_),
    filteredRandomBox_(ptf.filteredRandomBox_),
    U0_(ptf.U0_),
    computeVariant(ptf.computeVariant)
{}


Foam::turbulentDigitalFilterInletFvPatchVectorField::
turbulentDigitalFilterInletFvPatchVectorField
(
    const turbulentDigitalFilterInletFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    mapperPtr_(nullptr),
    variant_(ptf.variant_),
    isGaussian_(ptf.isGaussian_),
    isFixedSeed_(ptf.isFixedSeed_),
    isContinuous_(ptf.isContinuous_),
    isCorrectedFlowRate_(ptf.isCorrectedFlowRate_),
    interpolateR_(ptf.interpolateR_),
    interpolateUMean_(ptf.interpolateUMean_),
    isInsideMesh_(ptf.isInsideMesh_),
    isTaylorHypot_(ptf.isTaylorHypot_),
    mapMethod_(ptf.mapMethod_),
    curTimeIndex_(-1),
    tiny_(ptf.tiny_),
    patchNormalSpeed_(ptf.patchNormalSpeed_),
    modelConst_(ptf.modelConst_),
    perturb_(ptf.perturb_),
    initialFlowRate_(ptf.initialFlowRate_),
    rndGen_(ptf.rndGen_),
    planeDivisions_(ptf.planeDivisions_),
    invDelta_(ptf.invDelta_),
    indexPairs_(ptf.indexPairs_),
    R_(ptf.R_),
    LundWuSquires_(ptf.LundWuSquires_),
    UMean_(ptf.UMean_),
    Lbak_(ptf.Lbak_),
    L_(ptf.L_),
    const1FSM_(ptf.const1FSM_),
    const2FSM_(ptf.const2FSM_),
    constList1FSM_(ptf.constList1FSM_),
    constList2FSM_(ptf.constList2FSM_),
    lenRandomBox_(ptf.lenRandomBox_),
    randomBoxFactors2D_(ptf.randomBoxFactors2D_),
    randomBoxFactors3D_(ptf.randomBoxFactors3D_),
    iNextToLastPlane_(ptf.iNextToLastPlane_),
    randomBox_(ptf.randomBox_),
    filterCoeffs_(ptf.filterCoeffs_),
    filteredRandomBox_(ptf.filteredRandomBox_),
    U0_(ptf.U0_),
    computeVariant(ptf.computeVariant)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::turbulentDigitalFilterInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (curTimeIndex_ != db().time().timeIndex())
    {
        vectorField& U = *this;

        computeVariant(this, U);

        curTimeIndex_ = db().time().timeIndex();
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::turbulentDigitalFilterInletFvPatchVectorField::write
(
    Ostream& os
) const
{
    fixedValueFvPatchField<vector>::write(os),
    os.writeEntry("variant", variantNames[variant_]);
    os.writeEntry("planeDivisions", planeDivisions_);
    os.writeEntry("L", Lbak_);

    if (!interpolateR_)
    {
        R_.writeEntry("R", os);
    }

    if (!interpolateUMean_)
    {
        UMean_.writeEntry("UMean", os);
    }

    os.writeEntryIfDifferent<bool>("isGaussian", true, isGaussian_);
    os.writeEntryIfDifferent<bool>("isFixedSeed", true, isFixedSeed_);
    os.writeEntryIfDifferent<bool>("isContinuous", false, isContinuous_);
    os.writeEntryIfDifferent<bool>
        ("isCorrectedFlowRate", true, isCorrectedFlowRate_);
    os.writeEntryIfDifferent<bool>("isInsideMesh", false, isInsideMesh_);
    os.writeEntryIfDifferent<bool>("isTaylorHypot", true, isTaylorHypot_);

    if (!mapMethod_.empty())
    {
        os.writeEntryIfDifferent<word>
        (
            "mapMethod",
            "nearestCell",
            mapMethod_
        );
    }

    os.writeEntry("threshold", tiny_);
    os.writeEntry("patchNormalSpeed", patchNormalSpeed_);
    os.writeEntry("modelConst", modelConst_);
    os.writeEntry("perturb", perturb_);

    if (variant_ == variantType::FORWARD_STEPWISE)
    {
        os.writeEntry("const1FSM", const1FSM_);
        os.writeEntry("const2FSM", const2FSM_);
    }

    if (isContinuous_ && Pstream::master())
    {
        os.writeEntry("randomBox", randomBox_);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       turbulentDigitalFilterInletFvPatchVectorField
   );
}


// ************************************************************************* //
