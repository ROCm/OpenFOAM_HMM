/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::pointToPointPlanarInterpolation&
Foam::turbulentDigitalFilterInletFvPatchVectorField::patchMapper() const
{
    // Initialise interpolation (2D planar interpolation by triangulation)
    if (!mapperPtr_)
    {
        // Reread values and interpolate
        const fileName samplePointsFile
        (
            this->db().time().globalPath()
           /this->db().time().constant()
           /"boundaryData"
           /this->patch().name()
           /"points"
        );

        IOobject io
        (
            samplePointsFile,   // absolute path
            this->db().time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false,              // no need to register
            true                // is global object (currently not used)
        );

        // Read data
        const rawIOField<point> samplePoints(io, false);

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
Foam::turbulentDigitalFilterInletFvPatchVectorField::indexPairs()
{
    // Get patch normal direction into the domain
    const vector nf((-gAverage(patch().nf())).normalise());

    // Find the second local coordinate direction
    direction minCmpt = 0;
    scalar minMag = mag(nf[minCmpt]);
    for (direction cmpt = 1; cmpt < pTraits<vector>::nComponents; ++cmpt)
    {
        const scalar s = mag(nf[cmpt]);
        if (s < minMag)
        {
            minMag = s;
            minCmpt = cmpt;
        }
    }

    // Create the second local coordinate direction
    vector e2(Zero);
    e2[minCmpt] = 1;

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
    const Vector2D<label> n(n_.first(), n_.second());
    invDelta_[0] = n[0]/bb.span()[0];
    invDelta_[1] = n[1]/bb.span()[1];
    const point localMinPt(bb.min());

    // Compute virtual-actual patch index pairs
    List<Pair<label>> indexPairs(this->size(), Pair<label>(Zero, Zero));

    forAll(*this, facei)
    {
        const scalar& centre0 = localPos[facei][0];
        const scalar& centre1 = localPos[facei][1];

        // Virtual turbulence plane indices
        const label j = label((centre0 - localMinPt[0])*invDelta_[0]);
        const label k = label((centre1 - localMinPt[1])*invDelta_[1]);

        indexPairs[facei] = Pair<label>(facei, k*n[0] + j);
    }

    return indexPairs;
}


void Foam::turbulentDigitalFilterInletFvPatchVectorField::checkR() const
{
    const vectorField& faceCentres = this->patch().patch().faceCentres();

    forAll(R_, facei)
    {
        if (R_[facei].xx() <= 0)
        {
            FatalErrorInFunction
                << "Reynolds stress tensor component Rxx cannot be negative"
                << " or zero, where Rxx = " << R_[facei].xx()
                << " at the face centre = " << faceCentres[facei]
                << exit(FatalError);
        }

        if (R_[facei].yy() < 0 || R_[facei].zz() < 0)
        {
            FatalErrorInFunction
                << "Reynolds stress tensor components Ryy or Rzz cannot be"
                << " negative where Ryy = " << R_[facei].yy()
                << ", and Rzz = " << R_[facei].zz()
                << " at the face centre = " << faceCentres[facei]
                << exit(FatalError);
        }

        const scalar x0 = R_[facei].xx()*R_[facei].yy() - sqr(R_[facei].xy());

        if (x0 <= 0)
        {
            FatalErrorInFunction
                << "Reynolds stress tensor component group, Rxx*Ryy - Rxy^2"
                << " cannot be negative or zero"
                << " at the face centre = " << faceCentres[facei]
                << exit(FatalError);
        }

        const scalar x1 = R_[facei].zz() - sqr(R_[facei].xz())/R_[facei].xx();
        const scalar x2 = sqr(R_[facei].yz() - R_[facei].xy()*R_[facei].xz()
                         /(R_[facei].xx()*x0));
        const scalar x3 = x1 - x2;

        if (x3 < 0)
        {
            FatalErrorInFunction
                << "Reynolds stress tensor component group, "
                << "Rzz - Rxz^2/Rxx - (Ryz - Rxy*Rxz/(Rxx*(Rxx*Ryy - Rxy^2)))^2"
                << " cannot be negative at the face centre = "
                << faceCentres[facei]
                << exit(FatalError);
        }
    }

    Info<< "  # Reynolds stress tensor on patch is consistent #" << endl;
}


Foam::symmTensorField Foam::turbulentDigitalFilterInletFvPatchVectorField::
calcLund() const
{
    checkR();

    symmTensorField Lund(symmTensorField(R_.size()));

    forAll(R_, facei)
    {
        const symmTensor& R = R_[facei];
        symmTensor& lws = Lund[facei];

        // (KSJ:Eq. 5)
        lws.xx() = Foam::sqrt(R.xx());
        lws.xy() = R.xy()/lws.xx();
        lws.xz() = R.xz()/lws.xx();
        lws.yy() = Foam::sqrt(R.yy() - sqr(lws.xy()));
        lws.yz() = (R.yz() - lws.xy()*lws.xz())/lws.yy();
        lws.zz() = Foam::sqrt(R.zz() - sqr(lws.xz()) - sqr(lws.yz()));
    }

    return Lund;
}


Foam::scalar Foam::turbulentDigitalFilterInletFvPatchVectorField::
calcFlowRate() const
{
    const vector patchNormal((-gAverage(patch().nf())).normalise());
    return gSum((UMean_ & patchNormal)*patch().magSf());
}


Foam::tensor Foam::turbulentDigitalFilterInletFvPatchVectorField::
meterToCell
(
    const tensor& L
) const
{
    tensor Ls(L);

    // Convert streamwise integral length scales to integral
    // time scales by using Taylor's frozen turbulence hypothesis
    forAll(Ls.x(), i)
    {
        Ls[i] /= Ubulk_;
    }

    const scalar invDeltaT = 1.0/db().time().deltaTValue();

    //  (KSJ:Eq. 13)
    Ls.row(0, Ls.x()*invDeltaT);
    Ls.row(1, Ls.y()*invDelta_[0]);
    Ls.row(2, Ls.z()*invDelta_[1]);

    return Ls;
}


Foam::List<Foam::label> Foam::turbulentDigitalFilterInletFvPatchVectorField::
initBox() const
{
    label initValue = 0;
    label rangeModifier = 0;

    if (fsm_)
    {
        // Initialise with 1 since streamwise-dir possesses 1 cell in FSM
        initValue = 1;
        rangeModifier = pTraits<vector>::nComponents;
    }

    List<label> szBox(pTraits<tensor>::nComponents, initValue);
    Vector<label> szPlane
    (
        1,
        n_.first(),
        n_.second()
    );

    for
    (
        const label i
      : labelRange(rangeModifier, pTraits<tensor>::nComponents - rangeModifier)
    )
    {
        // Slicing index
        const label slicei = label(i/pTraits<vector>::nComponents);

        // Refer to "calcFilterCoeffs()"
        const label n = ceil(L_[i]);
        const label twiceN = 4*n;

        // Size of random-number sets
        szBox[i] = szPlane[slicei] + twiceN;
    }

    return szBox;
}


Foam::List<Foam::label> Foam::turbulentDigitalFilterInletFvPatchVectorField::
initFactors2D() const
{
    List<label> boxFactors2D(pTraits<vector>::nComponents);

    for (direction dir = 0; dir < pTraits<vector>::nComponents; ++dir)
    {
        boxFactors2D[dir] =
            szBox_[pTraits<vector>::nComponents + dir]
           *szBox_[pTraits<symmTensor>::nComponents + dir];
    }

    return boxFactors2D;
}


Foam::List<Foam::label> Foam::turbulentDigitalFilterInletFvPatchVectorField::
initFactors3D() const
{
    List<label> boxFactors3D(pTraits<vector>::nComponents);

    for (direction dir = 0; dir < pTraits<vector>::nComponents; ++dir)
    {
        boxFactors3D[dir] = boxFactors2D_[dir]*szBox_[dir];
    }

    return boxFactors3D;
}


Foam::List<Foam::label> Foam::turbulentDigitalFilterInletFvPatchVectorField::
initPlaneFactors() const
{
    List<label> planeFactors(pTraits<vector>::nComponents);

    for (direction dir = 0; dir < pTraits<vector>::nComponents; ++dir)
    {
        planeFactors[dir] = boxFactors2D_[dir]*(szBox_[dir] - 1);
    }

    return planeFactors;
}


Foam::List<Foam::List<Foam::scalar>>
Foam::turbulentDigitalFilterInletFvPatchVectorField::fillBox()
{
    List<List<scalar>> box(pTraits<vector>::nComponents, List<scalar>());

    // Initialise: Remaining convenience factors for (e1 e2 e3)
    for (direction dir = 0; dir < pTraits<vector>::nComponents; ++dir)
    {
        // Initialise: Random-box content with random-number sets
        box[dir] = randomSet<List<scalar>, scalar>(boxFactors3D_[dir]);
    }

    return box;
}


Foam::List<Foam::List<Foam::scalar>>
Foam::turbulentDigitalFilterInletFvPatchVectorField::calcFilterCoeffs() const
{
    List<List<scalar>> filterCoeffs
    (
        pTraits<tensor>::nComponents,
        List<scalar>(1, 1.0)
    );

    label initValue = 0;

    if (fsm_)
    {
        initValue = pTraits<vector>::nComponents;
    }

    for (direction dir = initValue; dir < pTraits<tensor>::nComponents; ++dir)
    {
        // The smallest filter width larger than length scale
        // (KSJ:'n' in Eq. 15)
        const label n  = ceil(L_[dir]);

        // Full filter-width (except mid-zero) according to the condition
        // (KSJ:Eq. 15 whereat N is minimum =2n)
        const label twiceN = 4*n;

        // Initialise filter-coeffs containers with full filter-width size
        // Extra elem is mid-zero within [-N, N]
        filterCoeffs[dir] = List<scalar>(twiceN + 1, Zero);

        // First element: -N within [-N, N]
        const scalar initElem = -2*scalar(n);

        // Container initialised with [-N, N] (KSJ:p. 658, item-b)
        std::iota
        (
            filterCoeffs[dir].begin(),
            filterCoeffs[dir].end(),
            initElem
        );

        // Compute filter-coeffs (KSJ:Eq. 14 (Gaussian))
        List<scalar> fTemp(filterCoeffs[dir]);
        scalar fSum = 0;
        const scalar nSqr = n*n;

        if (Gaussian_)
        {
            fTemp = sqr(exp(C1_*sqr(fTemp)/nSqr));
            fSum = Foam::sqrt(sum(fTemp));

            filterCoeffs[dir] =
                exp(C1_*sqr(filterCoeffs[dir])/nSqr)/fSum;
        }
        else
        {
            fTemp = sqr(exp(C1_*mag(fTemp)/n));
            fSum = Foam::sqrt(sum(fTemp));

            filterCoeffs[dir] = exp(C1_*mag(filterCoeffs[dir])/n)/fSum;
        }
    }

    return filterCoeffs;
}


void Foam::turbulentDigitalFilterInletFvPatchVectorField::shiftRefill()
{
    forAll(box_, dir)
    {
        List<scalar>& r = box_[dir];

        // Shift forward from the back to the front / First Out
        inplaceRotateList(r, boxFactors2D_[dir]);

        // Refill the back with a new random-number set / First In
        for (label i = 0; i < boxFactors2D_[dir]; ++i)
        {
            r[i] = rndGen_.GaussNormal<scalar>();
        }
    }
}


void Foam::turbulentDigitalFilterInletFvPatchVectorField::mapFilteredBox
(
    vectorField& U
)
{
    for (const auto& x : indexPairs_)
    {
        const label xf = x.first();
        const label xs = x.second();
        U[xf][0] = filteredBox_[0][xs];
        U[xf][1] = filteredBox_[1][xs];
        U[xf][2] = filteredBox_[2][xs];
    }
}


void Foam::turbulentDigitalFilterInletFvPatchVectorField::onePointCorrs
(
    vectorField& U
) const
{
    forAll(Lund_, facei)
    {
        vector& Us = U[facei];
        const symmTensor& lws = Lund_[facei];

        // (KSJ:p. 658, item-e)
        Us.z() = Us.x()*lws.xz() + Us.y()*lws.yz() + Us.z()*lws.zz();
        Us.y() = Us.x()*lws.xy() + Us.y()*lws.yy();
        Us.x() = Us.x()*lws.xx();
    }
}


void Foam::turbulentDigitalFilterInletFvPatchVectorField::twoPointCorrs()
{
    for (direction dir = 0; dir < pTraits<vector>::nComponents; ++dir)
    {
        List<scalar>& in = box_[dir];
        List<scalar>& out = filteredBox_[dir];
        const List<scalar>& filter1 = filterCoeffs_[dir];
        const List<scalar>& filter2 = filterCoeffs_[3 + dir];
        const List<scalar>& filter3 = filterCoeffs_[6 + dir];

        const label sz1 = szBox_[dir];
        const label sz2 = szBox_[3 + dir];
        const label sz3 = szBox_[6 + dir];
        const label szfilter1 = filterCoeffs_[dir].size();
        const label szfilter2 = filterCoeffs_[3 + dir].size();
        const label szfilter3 = filterCoeffs_[6 + dir].size();
        const label sz23 = boxFactors2D_[dir];
        const label sz123 = boxFactors3D_[dir];
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


Foam::List<Foam::scalar> Foam::turbulentDigitalFilterInletFvPatchVectorField::
calcCoeffs1FSM() const
{
    List<scalar> coeffs1FSM(pTraits<vector>::nComponents);

    forAll(L_.x(), i)
    {
        coeffs1FSM[i] = exp(C1FSM_/L_[i]);
    }

    return coeffs1FSM;
}


Foam::List<Foam::scalar> Foam::turbulentDigitalFilterInletFvPatchVectorField::
calcCoeffs2FSM() const
{
    List<scalar> coeffs2FSM(pTraits<vector>::nComponents);

    forAll(L_.x(), i)
    {
        coeffs2FSM[i] = Foam::sqrt(1.0 - exp(C2FSM_/L_[i]));
    }

    return coeffs2FSM;
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
    fsm_(false),
    Gaussian_(true),
    fixSeed_(true),
    continuous_(false),
    correctFlowRate_(true),
    interpR_(false),
    interpUMean_(false),
    mapMethod_("nearestCell"),
    curTimeIndex_(-1),
    Ubulk_(0.0),
    C1_(-0.5*constant::mathematical::pi),
    perturb_(1e-5),
    flowRate_(1.0),
    rndGen_(1234),
    n_(0, 0),
    invDelta_(Zero),
    indexPairs_(Zero),
    R_(Zero),
    Lund_(Zero),
    UMean_(Zero),
    Lbak_(Zero),
    L_(Zero),
    C1FSM_(Zero),
    C2FSM_(Zero),
    coeffs1FSM_(Zero),
    coeffs2FSM_(Zero),
    szBox_(Zero),
    boxFactors2D_(Zero),
    boxFactors3D_(Zero),
    iNextToLastPlane_(Zero),
    box_(Zero),
    filterCoeffs_(Zero),
    filteredBox_(Zero),
    U0_(Zero)
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
    fsm_(ptf.fsm_),
    Gaussian_(ptf.Gaussian_),
    fixSeed_(ptf.fixSeed_),
    continuous_(ptf.continuous_),
    correctFlowRate_(ptf.correctFlowRate_),
    interpR_(ptf.interpR_),
    interpUMean_(ptf.interpUMean_),
    mapMethod_(ptf.mapMethod_),
    curTimeIndex_(-1),
    Ubulk_(ptf.Ubulk_),
    C1_(ptf.C1_),
    perturb_(ptf.perturb_),
    flowRate_(ptf.flowRate_),
    rndGen_(ptf.rndGen_),
    n_(ptf.n_),
    invDelta_(ptf.invDelta_),
    indexPairs_(ptf.indexPairs_),
    R_(ptf.R_),
    Lund_(ptf.Lund_),
    UMean_(ptf.UMean_),
    Lbak_(ptf.Lbak_),
    L_(ptf.L_),
    C1FSM_(ptf.C1FSM_),
    C2FSM_(ptf.C2FSM_),
    coeffs1FSM_(ptf.coeffs1FSM_),
    coeffs2FSM_(ptf.coeffs2FSM_),
    szBox_(ptf.szBox_),
    boxFactors2D_(ptf.boxFactors2D_),
    boxFactors3D_(ptf.boxFactors3D_),
    iNextToLastPlane_(ptf.iNextToLastPlane_),
    box_(ptf.box_),
    filterCoeffs_(ptf.filterCoeffs_),
    filteredBox_(ptf.filteredBox_),
    U0_(ptf.U0_)
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
    fsm_(dict.getOrDefault("fsm", false)),
    Gaussian_(fsm_ ? false : dict.getOrDefault("Gaussian", true)),
    fixSeed_(dict.getOrDefault("fixSeed", true)),
    continuous_(dict.getOrDefault("continuous", false)),
    correctFlowRate_(dict.getOrDefault("correctFlowRate", true)),
    interpR_(false),
    interpUMean_(false),
    mapMethod_(dict.getOrDefault<word>("mapMethod", "nearestCell")),
    curTimeIndex_(-1),
    Ubulk_(dict.getCheck<scalar>("Ubulk", scalarMinMax::ge(ROOTSMALL))),
    C1_
    (
        dict.getOrDefault<scalar>("C1", -0.5*constant::mathematical::pi)
    ),
    perturb_
    (
        dict.getCheckOrDefault<scalar>("perturb", 1e-5, scalarMinMax::ge(SMALL))
    ),
    flowRate_(1.0),
    rndGen_(fixSeed_ ? 1234 : time(0)),
    n_
    (
        dict.getCheck<Tuple2<label, label>>
        (
            "n",
            [&](const Tuple2<label, label>& x)
            {
                return min(x.first(), x.second()) > 0 ? true : false;
            }
        )
    ),
    invDelta_(Zero),
    indexPairs_(indexPairs()),
    R_(interpolateOrRead<symmTensor>("R", dict, interpR_)),
    Lund_(calcLund()),
    UMean_(interpolateOrRead<vector>("UMean", dict, interpUMean_)),
    Lbak_
    (
        dict.getCheck<tensor>
        (
            "L",
            [&](const tensor& x){return cmptMin(x) > ROOTSMALL ? true : false;}
        )
    ),
    L_(meterToCell(Lbak_)),
    C1FSM_
    (
        fsm_
      ? dict.getOrDefault<scalar>
        (
            "C1FSM",
           -0.25*constant::mathematical::pi
        )
      : 0.0
    ),
    C2FSM_
    (
        fsm_
      ? dict.getOrDefault<scalar>
        (
            "C2FSM",
            -0.5*constant::mathematical::pi
        )
      : 0.0
    ),
    coeffs1FSM_(fsm_ ? calcCoeffs1FSM() : List<scalar>()),
    coeffs2FSM_(fsm_ ? calcCoeffs2FSM() : List<scalar>()),
    szBox_(initBox()),
    boxFactors2D_(initFactors2D()),
    boxFactors3D_(initFactors3D()),
    iNextToLastPlane_(initPlaneFactors()),
    box_
    (
        (continuous_ && Pstream::master())
      ? dict.getOrDefault<List<List<scalar>>>
        (
            "box",
            fillBox() // First time-step
        )
      :
        fillBox()
    ),
    filterCoeffs_(calcFilterCoeffs()),
    filteredBox_
    (
        pTraits<vector>::nComponents,
        List<scalar>(n_.first()*n_.second(), Zero)
    ),
    U0_
    (
        fsm_
      ? randomSet<vectorField, vector>(patch().size())
      : vectorField()
    )
{
    if (correctFlowRate_)
    {
        flowRate_ = calcFlowRate();
    }

    // Check if varying or fixed time-step computation
    if (db().time().isAdjustTimeStep())
    {
        WarningInFunction
            << "  # Varying time-step computations are not fully supported #"
            << endl;
    }

    Info<< "  # turbulentDigitalFilterInlet initialisation completed #" << endl;
}


Foam::turbulentDigitalFilterInletFvPatchVectorField::
turbulentDigitalFilterInletFvPatchVectorField
(
    const turbulentDigitalFilterInletFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    mapperPtr_(nullptr),
    fsm_(ptf.fsm_),
    Gaussian_(ptf.Gaussian_),
    fixSeed_(ptf.fixSeed_),
    continuous_(ptf.continuous_),
    correctFlowRate_(ptf.correctFlowRate_),
    interpR_(ptf.interpR_),
    interpUMean_(ptf.interpUMean_),
    mapMethod_(ptf.mapMethod_),
    curTimeIndex_(ptf.curTimeIndex_),
    Ubulk_(ptf.Ubulk_),
    C1_(ptf.C1_),
    perturb_(ptf.perturb_),
    flowRate_(ptf.flowRate_),
    rndGen_(ptf.rndGen_),
    n_(ptf.n_),
    invDelta_(ptf.invDelta_),
    indexPairs_(ptf.indexPairs_),
    R_(ptf.R_),
    Lund_(ptf.Lund_),
    UMean_(ptf.UMean_),
    Lbak_(ptf.Lbak_),
    L_(ptf.L_),
    C1FSM_(ptf.C1FSM_),
    C2FSM_(ptf.C2FSM_),
    coeffs1FSM_(ptf.coeffs1FSM_),
    coeffs2FSM_(ptf.coeffs2FSM_),
    szBox_(ptf.szBox_),
    boxFactors2D_(ptf.boxFactors2D_),
    boxFactors3D_(ptf.boxFactors3D_),
    iNextToLastPlane_(ptf.iNextToLastPlane_),
    box_(ptf.box_),
    filterCoeffs_(ptf.filterCoeffs_),
    filteredBox_(ptf.filteredBox_),
    U0_(ptf.U0_)
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
    fsm_(ptf.fsm_),
    Gaussian_(ptf.Gaussian_),
    fixSeed_(ptf.fixSeed_),
    continuous_(ptf.continuous_),
    correctFlowRate_(ptf.correctFlowRate_),
    interpR_(ptf.interpR_),
    interpUMean_(ptf.interpUMean_),
    mapMethod_(ptf.mapMethod_),
    curTimeIndex_(-1),
    Ubulk_(ptf.Ubulk_),
    C1_(ptf.C1_),
    perturb_(ptf.perturb_),
    flowRate_(ptf.flowRate_),
    rndGen_(ptf.rndGen_),
    n_(ptf.n_),
    invDelta_(ptf.invDelta_),
    indexPairs_(ptf.indexPairs_),
    R_(ptf.R_),
    Lund_(ptf.Lund_),
    UMean_(ptf.UMean_),
    Lbak_(ptf.Lbak_),
    L_(ptf.L_),
    C1FSM_(ptf.C1FSM_),
    C2FSM_(ptf.C2FSM_),
    coeffs1FSM_(ptf.coeffs1FSM_),
    coeffs2FSM_(ptf.coeffs2FSM_),
    szBox_(ptf.szBox_),
    boxFactors2D_(ptf.boxFactors2D_),
    boxFactors3D_(ptf.boxFactors3D_),
    iNextToLastPlane_(ptf.iNextToLastPlane_),
    box_(ptf.box_),
    filterCoeffs_(ptf.filterCoeffs_),
    filteredBox_(ptf.filteredBox_),
    U0_(ptf.U0_)
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

        if (Pstream::master())
        {
            twoPointCorrs();
            shiftRefill();
        }

        Pstream::scatter(filteredBox_);

        mapFilteredBox(U);

        if (fsm_)
        {
            //(XC:Eq. 14)
            for (direction dir = 0; dir < pTraits<vector>::nComponents; ++dir)
            {
                U.component(dir) =
                    U0_.component(dir)*coeffs1FSM_[dir]
                  + U.component(dir)*coeffs2FSM_[dir];
            }

            U0_ = U;
        }

        onePointCorrs(U);

        U += UMean_;

        if (correctFlowRate_)
        {
            // (KCX:Eq. 8)
            U *= (flowRate_/gSum(U & -patch().Sf()));
        }

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
    os.writeEntry("n", n_);
    os.writeEntry("L", Lbak_);

    if (!interpR_)
    {
        R_.writeEntry("R", os);
    }

    if (!interpUMean_)
    {
        UMean_.writeEntry("UMean", os);
    }

    os.writeEntryIfDifferent<bool>("fsm", false, fsm_);
    os.writeEntryIfDifferent<bool>("Gaussian", true, Gaussian_);
    os.writeEntryIfDifferent<bool>("fixSeed", true, fixSeed_);
    os.writeEntryIfDifferent<bool>("continuous", false, continuous_);
    os.writeEntryIfDifferent<bool>("correctFlowRate", true, correctFlowRate_);

    if (!mapMethod_.empty())
    {
        os.writeEntryIfDifferent<word>
        (
            "mapMethod",
            "nearestCell",
            mapMethod_
        );
    }

    os.writeEntry("Ubulk", Ubulk_);
    os.writeEntry("C1", C1_);
    os.writeEntry("perturb", perturb_);

    if (fsm_)
    {
        os.writeEntry("C1FSM", C1FSM_);
        os.writeEntry("C2FSM", C2FSM_);
    }

    if (continuous_ && Pstream::master())
    {
        os.writeEntry("box", box_);
    }
}


void Foam::turbulentDigitalFilterInletFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<vector>::autoMap(m);
}


void Foam::turbulentDigitalFilterInletFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<vector>::rmap(ptf, addr);
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
