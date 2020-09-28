/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "STDMD.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(STDMD, 0);
    addToRunTimeSelectionTable(functionObject, STDMD, dictionary);
}
}


const Foam::Enum
<
    Foam::functionObjects::STDMD::modeSorterType
>
Foam::functionObjects::STDMD::modeSorterTypeNames
({
    { modeSorterType::KIEWAT , "kiewat" },
    { modeSorterType::KOU_ZHANG , "kouZhang" },
    { modeSorterType::FIRST_SNAPSHOT, "firstSnap" }
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::functionObjects::STDMD::parnorm
(
    const RMatrix& colVector
) const
{
    #ifdef FULLDEBUG
    // Check if the given RectangularMatrix is effectively a column vector
    if (colVector.n() != 1)
    {
        FatalErrorInFunction
            << "  # Input matrix is not a column vector. #"
            << abort(FatalError);
    }
    #endif
    const bool noSqrt = true;
    scalar result(colVector.columnNorm(0, noSqrt));
    reduce(result, sumOp<scalar>());

    // Heuristic addition to avoid very small or zero norm
    return max(SMALL, Foam::sqrt(result));
}


void Foam::functionObjects::STDMD::snapshot()
{
    bool processed = false;
    processed = processed || getSnapshotType<scalar>();
    processed = processed || getSnapshotType<vector>();
    processed = processed || getSnapshotType<sphericalTensor>();
    processed = processed || getSnapshotType<symmTensor>();
    processed = processed || getSnapshotType<tensor>();

    if (!processed)
    {
        FatalErrorInFunction
            << "  # Unknown type of input field during snapshot loading = "
            << fieldName_ << " #" << nl
            << "  # Do you execute required functionObjects "
            << "before executing STDMD, e.g. mapFields?"
            << abort(FatalError);
    }
}


void Foam::functionObjects::STDMD::init()
{
    bool processed = false;
    processed = processed || getComps<scalar>();
    processed = processed || getComps<vector>();
    processed = processed || getComps<sphericalTensor>();
    processed = processed || getComps<symmTensor>();
    processed = processed || getComps<tensor>();

    if (!processed)
    {
        FatalErrorInFunction
            << "  # Unknown type of input field during initialisation = "
            << fieldName_ << " #" << nl
            << "  # Do you execute required functionObjects "
            << "before executing STDMD, e.g. mapFields?"
            << abort(FatalError);
    }

    nSnap_ = nComps_*mesh_.nCells();

    if (nSnap_ <= 0)
    {
        FatalErrorInFunction
            << "  # Zero-size input field = " << fieldName_ << " #"
            << abort(FatalError);
    }

    currIndex_ = 0;
    zNorm_ = 0;
    ezNorm_ = 0;

    z_ = RMatrix(2*nSnap_, 1, Zero);
    ez_ = z_;
    X1_ = RMatrix(nSnap_, 1, Zero);
    Qz_ = z_;
    Gz_ = SMatrix(1);

    RxInv_.clear();
    Ap_.clear();
    oEVecs_.clear();
    oEVals_.clear();
    oAmps_.clear();
    oFreqs_.clear();
    iFreqs_.clear();
    oMags_.clear();
    iMags_.clear();

    initialised_ = true;
}


void Foam::functionObjects::STDMD::initBasis()
{
    std::copy(z_.cbegin(), z_.cbegin() + nSnap_, X1_.begin());

    zNorm_ = parnorm(z_);
    Qz_ = z_/zNorm_;
    Gz_(0,0) = sqr(zNorm_);
}


void Foam::functionObjects::STDMD::GramSchmidt()
{
    ez_ = z_;
    RMatrix dz(Qz_.n(), 1, Zero);

    for (label i = 0; i < nGramSchmidt_; ++i)
    {
        dz = Qz_ & ez_;
        reduce(dz, sumOp<RMatrix>());
        ez_ -= Qz_*dz;
    }
}


void Foam::functionObjects::STDMD::expandBasis()
{
    Log<< tab << "# " << name() << ":"
        << " Expanding orthonormal basis for field = " << fieldName_ << " #"
        << endl;

    // Stack a column (ez_/ezNorm_) to Qz_
    Qz_.resize(Qz_.m(), Qz_.n() + 1);
    Qz_.subColumn(Qz_.n() - 1) = ez_/ezNorm_;

    // Stack a row (Zero) and column (Zero) to Gz_
    Gz_.resize(Gz_.m() + 1);
}


void Foam::functionObjects::STDMD::updateGz()
{
    RMatrix zTilde(Qz_ & z_);
    reduce(zTilde, sumOp<RMatrix>());

    const SMatrix zTildes(zTilde^zTilde);

    Gz_ += zTildes;
}


void Foam::functionObjects::STDMD::compressBasis()
{
    Log<< tab << "# " << name() << ":"
        << " Compressing orthonormal basis for field = " << fieldName_ << " #"
        << endl;

    RMatrix qz;

    if (Pstream::master())
    {
        const bool symmetric = true;
        const EigenMatrix<scalar> EM(Gz_, symmetric);
        const SquareMatrix<scalar>& EVecs = EM.EVecs();
        DiagonalMatrix<scalar> EVals(EM.EValsRe());

        if (testEigen_)
        {
            testEigenvalues(Gz_, EVals);
            testEigenvectors(Gz_, EVals, EVecs);
        }

        // Sort eigenvalues in descending order, and track indices
        const auto descend = [&](scalar a, scalar b){ return a > b; };
        const List<label> permut(EVals.sortPermutation(descend));
        EVals.applyPermutation(permut);
        EVals.resize(EVals.size() - 1);

        // Update Gz_
        Gz_ = SMatrix(maxRank_, Zero);
        Gz_.diag(EVals);

        qz.resize(Qz_.n(), maxRank_);
        for (label i = 0; i < maxRank_; ++i)
        {
            qz.subColumn(i) = EVecs.subColumn(permut[i]);
        }
    }
    Pstream::scatter(Gz_);
    Pstream::scatter(qz);

    // Update Qz_
    Qz_ = Qz_*qz;
}


void Foam::functionObjects::STDMD::calcAp()
{
    Log<< tab << "# " << name() << ": Computing Ap matrix #" << endl;

    Log<< tab << "# " << name() << ": Computing local Rx #" << endl;
    RMatrix Rx;
    {
        QRMatrix<RMatrix> QRM
        (
            QzUH_,
            QRMatrix<RMatrix>::outputTypes::REDUCED_R,
            QRMatrix<RMatrix>::colPivoting::FALSE
        );
        Rx = QRM.R();
    }
    Rx.round();

    RMatrix A1; // Convenience objects A1, A2, A3 to compute Ap
    if (Pstream::parRun())
    {
        // Parallel direct tall-skinny QR decomposition
        // (BGD:Fig. 5) & (DGHL:Fig. 2)
        // Tests revealed that the distribution of Qz_ does not affect
        // the final outcome of TSQR decomposition up to sign

        Log<< tab << "# " << name() << ": Gathering all local Rx #" << endl;
        List<RMatrix> RxList(Pstream::nProcs());
        RxList[Pstream::myProcNo()] = Rx;
        Pstream::gatherList(RxList);
        Pstream::scatterList(RxList);

        // Create a global object Rx
        if (Pstream::master())
        {
            // Determine the row size of the global Rx
            label mRowSum = 0;
            forAll(RxList, i)
            {
                mRowSum += RxList[i].m();
            }

            Rx.resize(mRowSum, Rx.n());

            Log<< tab << "# " << name()
                << ": Populating the global Rx #" << endl;
            label m = 0;
            for (const int i : Pstream::allProcs())
            {
                const label mRows = RxList[i].m();

                Rx.subMatrix(m, 0, mRows) = RxList[i];

                m += mRows;
            }

            Log<< tab << "# " << name()
                << ": Computing the parallel-direct tall-skinny QR decomp. #"
                << endl;
            QRMatrix<RMatrix> QRM
            (
                Rx,
                QRMatrix<RMatrix>::outputTypes::REDUCED_R,
                QRMatrix<RMatrix>::storeMethods::IN_PLACE,
                QRMatrix<RMatrix>::colPivoting::FALSE
            );
            Rx.round();

            Log<< tab << "# " << name()
                << ": Computing Moore-Penrose pseudo-inverse of Rx #"
                << endl;
            RxInv_ = pinv(Rx);

            Log<< tab << "# " << name() << ": Computing Gx #" << endl;
            const RMatrix Gx(Rx*(Gz_^Rx));

            Log<< tab << "# " << name()
                << ": Computing Moore-Penrose pseudo-inverse of Gx #"
                << endl;
            const RMatrix GxInv(pinv(Gx));

            Log<< tab << "# " << name() << ": Computing A1 #" << endl;
            A1 = RxInv_*GxInv;
        }
        Pstream::scatter(RxInv_);
        Pstream::scatter(A1);

        Log<< tab << "# " << name() << ": Computing A2 #" << endl;
        SMatrix A2(QzUH_ & QzUH_);
        reduce(A2, sumOp<SMatrix>());

        Log<< tab << "# " << name() << ": Computing A3 #" << endl;
        SMatrix A3(QzUH_ & QzLH_);
        reduce(A3, sumOp<SMatrix>());

        Log<< tab << "# " << name() << ": Computing Ap #" << endl;
        // by optimized matrix chain multiplication
        // obtained by dynamic programming memoisation
        Ap_ = SMatrix(RxInv_ & ((A3*(Gz_*A2))*A1));
    }
    else
    {
        {
            Log<< tab << "# " << name()
                << ": Computing Moore-Penrose pseudo-inverse of Rx #"
                << endl;
            RxInv_ = pinv(Rx);

            Log<< tab << "# " << name() << ": Computing Gx #" << endl;
            const RMatrix Gx(Rx*(Gz_^Rx));

            Log<< tab << "# " << name()
                << ": Computing Moore-Penrose pseudo-inverse of Gx #"
                << endl;
            const RMatrix GxInv(pinv(Gx));

            Log<< tab << "# " << name() << ": Computing A1 #" << endl;
            A1 = RxInv_*GxInv;
        }

        Log<< tab << "# " << name() << ": Computing A2 #" << endl;
        SMatrix A2(QzUH_ & QzUH_);

        Log<< tab << "# " << name() << ": Computing A3 #" << endl;
        SMatrix A3(QzUH_ & QzLH_);

        Log<< tab << "# " << name() << ": Computing Ap #" << endl;
        // by optimized matrix chain multiplication
        // obtained by dynamic programming memoisation
        Ap_ = SMatrix(RxInv_ & ((A3*(Gz_*A2))*A1));
    }
}


void Foam::functionObjects::STDMD::calcEigen()
{
    Log<< tab << "# " << name() << ": Computing eigendecomposition #" << endl;

    // Compute eigenvalues, and clip eigenvalues with values < minMagEVal_
    if (Pstream::master())
    {
        // Replace Ap by eigenvalues (in-place eigendecomposition)
        const EigenMatrix<scalar> EM(Ap_);
        const DiagonalMatrix<scalar>& evalsRe = EM.EValsRe();
        const DiagonalMatrix<scalar>& evalsIm = EM.EValsIm();

        oEVals_.resize(evalsRe.size());
        oEVecs_ = RectangularMatrix<complex>(EM.complexEVecs());

        // Convert scalar eigenvalue pairs to complex eigenvalues
        label i = 0;
        for (auto& eval : oEVals_)
        {
            eval = complex(evalsRe[i], evalsIm[i]);
            ++i;
        }

        if (testEigen_)
        {
            testEigenvalues(Ap_, evalsRe);
            testEigenvectors(Ap_, oEVals_, oEVecs_);
        }
    }
}


bool Foam::functionObjects::STDMD::close
(
    const scalar s1,
    const scalar s2,
    const scalar absTol,
    const scalar relTol
) const
{
    if (s1 == s2) return true;

    return (mag(s1 - s2) <= max(relTol*max(mag(s1), mag(s2)), absTol));
}


void Foam::functionObjects::STDMD::testEigenvalues
(
    const SquareMatrix<scalar>& A,
    const DiagonalMatrix<scalar>& EValsRe
) const
{
    const scalar trace = A.trace();
    // Imaginary part of complex conjugates cancel each other
    const scalar EValsSum = sum(EValsRe);

    const bool clse = close(EValsSum - trace, 0, absTol_, relTol_);

    if (clse)
    {
        Info<< tab
            << "  ## PASS: Eigenvalues ##" << endl;
    }
    else
    {
        Info<< tab
            << "  ## INCONCLUSIVE: Eigenvalues ##" << nl
            << "  # sum(eigenvalues) = " << EValsSum << nl
            << "  # trace(A) = " << trace << nl
            << "  # sum(eigenvalues) - trace(A) ~ 0 ?= "
            << EValsSum - trace << nl
            << "  #######################"
            << endl;

        if (dumpEigen_)
        {
            Info<< tab
                << "  ## Operands ##" << nl
                << "  # eigenvalues:" << nl << EValsRe << nl
                << "  # input matrix A:" << nl << A << nl
                << "  ##############"
                << endl;
        }
    }
}


void Foam::functionObjects::STDMD::testEigenvectors
(
    const SquareMatrix<scalar>& A,
    const DiagonalMatrix<scalar>& EValsRe,
    const SquareMatrix<scalar>& EVecs
) const
{
    unsigned nInconclusive = 0;

    for (label i = 0; i < A.n(); ++i)
    {
        const RectangularMatrix<scalar> EVec(EVecs.subColumn(i));
        const scalar EVal = EValsRe[i];
        const RectangularMatrix<scalar> leftSide(A*EVec);
        const RectangularMatrix<scalar> rightSide(EVal*EVec);
        const scalar rslt = (leftSide - rightSide).norm();

        const bool clse = close(rslt, 0, absTol_, relTol_);

        if (!clse)
        {
            Info<< tab
                << "  ## INCONCLUSIVE: Eigenvector ##" << nl
                << "  # (A & EVec - EVal*EVec).norm() ~ 0 ?= " << rslt << nl
                << "  ##################################"
                << endl;

            if (dumpEigen_)
            {
                Info<< tab
                    << "  ## Operands ##" << nl
                    << "  # eigenvalue:" << nl << EVal << nl
                    << "  # input matrix A:" << nl << A << nl
                    << "  # eigenvector:" << nl << EVec << nl
                    << "  # (A & EVec):" << nl << leftSide << nl
                    << "  # (EVal*EVec):" << nl << rightSide << nl
                    << "  ##############"
                    << endl;
            }

            ++nInconclusive;
        }
    }

    if (!nInconclusive)
    {
        Info<< tab
            << "  ## PASS: Eigenvectors ##" << endl;
    }
}


void Foam::functionObjects::STDMD::testEigenvectors
(
    const SquareMatrix<scalar>& A,
    const List<complex>& EVals,
    const RectangularMatrix<complex>& EVecs
) const
{
    SquareMatrix<complex> B(A.m());
    auto convertToComplex = [&](const scalar& val) { return complex(val); };
    std::transform
    (
        A.cbegin(),
        A.cend(),
        B.begin(),
        convertToComplex
    );

    unsigned nInconclusive = 0;

    for (label i = 0; i < B.n(); ++i)
    {
        const RectangularMatrix<complex> EVec(EVecs.subColumn(i));
        const complex EVal = EVals[i];
        const RectangularMatrix<complex> leftSide(B*EVec);
        const RectangularMatrix<complex> rightSide(EVal*EVec);
        const scalar rslt = (leftSide - rightSide).norm();

        const bool clse = close(rslt, 0, absTol_, relTol_);

        if (!clse)
        {
            Info<< tab
                << "  ## INCONCLUSIVE: Eigenvector ##" << nl
                << "  # (A & EVec - EVal*EVec).norm() ~ 0 ?= " << rslt << nl
                << "  ##################################"
                << endl;

            if (dumpEigen_)
            {
                Info<< tab
                    << "  ## Operands ##" << nl
                    << "  # eigenvalue:" << nl << EVal << nl
                    << "  # input matrix A:" << nl << A << nl
                    << "  # eigenvector:" << nl << EVec << nl
                    << "  # (A & EVec):" << nl << leftSide << nl
                    << "  # (EVal*EVec):" << nl << rightSide << nl
                    << "  ##############"
                    << endl;
            }

            ++nInconclusive;
        }
    }

    if (!nInconclusive)
    {
        Info<< tab
            << "  ## PASS: Eigenvectors ##" << endl;
    }
}


void Foam::functionObjects::STDMD::filterEVals()
{
    Log<< tab << "# " << name() << ": Filtering eigenvalues #" << endl;

    if (Pstream::master())
    {
        List<complex> cpEVals(oEVals_.size());
        auto it =
            std::copy_if
            (
                oEVals_.cbegin(),
                oEVals_.cend(),
                cpEVals.begin(),
                [&](const complex& x){ return minMagEVal_ < mag(x); }
            );
        cpEVals.resize(std::distance(cpEVals.begin(), it));

        if (cpEVals.size() == 0)
        {
            WarningInFunction
                << "No eigenvalue with mag(eigenvalue) larger than "
                << "minMagEVal_ = " << minMagEVal_ << " was found."
                << endl;
        }
        else
        {
            oEVals_ = cpEVals;
        }
    }
    Pstream::scatter(oEVals_);
    Pstream::scatter(oEVecs_);
}


void Foam::functionObjects::STDMD::calcFreqs()
{
    Log<< tab << "# " << name() << ": Computing frequencies #" << endl;

    if (Pstream::master())
    {
        oFreqs_.resize(oEVals_.size());

        // Frequency equation (K:Eq. 81)
        auto fEq =
            [&](const complex& EVal)
            {
                return Foam::log(max(EVal, complex(SMALL))).imag()/(twoPi*dt_);
            };

        // Compute frequencies
        std::transform(oEVals_.cbegin(), oEVals_.cend(), oFreqs_.begin(), fEq);
    }
}


void Foam::functionObjects::STDMD::calcFreqI()
{
    Log<< tab << "# " << name() << ": Computing frequency indices #" << endl;

    if (Pstream::master())
    {
        // Initialise iterator by the first search
        auto margin = [&](const scalar& x){ return (fMin_ <= x && x < fMax_); };

        auto it = std::find_if(oFreqs_.cbegin(), oFreqs_.cend(), margin);

        while (it != oFreqs_.end())
        {
            iFreqs_.append(std::distance(oFreqs_.cbegin(), it));
            it = std::find_if(std::next(it), oFreqs_.cend(), margin);
        }
    }
    Pstream::scatter(oFreqs_);
    Pstream::scatter(iFreqs_);
}


void Foam::functionObjects::STDMD::calcAmps()
{
    Log<< tab << "# " << name() << ": Computing amplitudes #" << endl;

    RMatrix temp((RxInv_.T()^QzUH_)*X1_);
    reduce(temp, sumOp<RMatrix>());

    if (Pstream::master())
    {
        oAmps_.resize(temp.m());
        const RCMatrix pinvEVecs(pinv(oEVecs_));

        // oAmps_ = pinvEVecs*temp;
        for (label i = 0; i < oAmps_.size(); ++i)
        {
            for (label k = 0; k < temp.m(); ++k)
            {
                oAmps_[i] += pinvEVecs(i, k)*temp(k, 0);
            }
        }
    }
    Pstream::scatter(oAmps_);
}


void Foam::functionObjects::STDMD::calcMags()
{
    Log<< tab << "# " << name() << ": Computing magnitudes #" << endl;

    if (Pstream::master())
    {
        oMags_.resize(oAmps_.size());

        Log<< tab << "# " << name() << ": Sorting modes with ";

        switch (modeSorter_)
        {
            case modeSorterType::FIRST_SNAPSHOT:
            {
                Log<< "method of first snapshot #" << endl;

                std::transform
                (
                    oAmps_.cbegin(),
                    oAmps_.cend(),
                    oMags_.begin(),
                    [&](const complex& val){ return mag(val); }
                );
                break;
            }

            case modeSorterType::KIEWAT:
            {
                Log<< "modified weighted amplitude scaling method #" << endl;

                //  Eigendecomposition returns eigenvectors with
                //  the unity norm, hence modeNorm = 1
                const scalar modeNorm = 1;
                const scalar pr = 1;
                List<scalar> w(currIndex_);
                std::iota(w.begin(), w.end(), 1);
                w = sin(twoPi/currIndex_*(w - 1 - 0.25*currIndex_))*pr + pr;

                label i = 0;
                for (const auto& amp : oAmps_)
                {
                    oMags_[i] = sorter(w, amp, oEVals_[i], modeNorm);
                    ++i;
                }
                break;
            }

            case modeSorterType::KOU_ZHANG:
            {
                Log<< "weighted amplitude scaling method #" << endl;

                const scalar modeNorm = 1;
                const List<scalar> w(currIndex_, 1.0);
                label i = 0;
                for (const auto& amp : oAmps_)
                {
                    oMags_[i] = sorter(w, amp, oEVals_[i], modeNorm);
                    ++i;
                }
                break;
            }

            default:
                break;
        }
    }
}


void Foam::functionObjects::STDMD::calcMagI()
{
    Log<< tab << "# " << name() << ": Computing magnitude indices #" << endl;

    if (Pstream::master())
    {
        iMags_ = iFreqs_;

        auto descend =
            [&](const label i1, const label i2)
            {
                return !(oMags_[i1] < oMags_[i2]);
            };

        std::sort(iMags_.begin(), iMags_.end(), descend);
    }
    Pstream::scatter(oMags_);
    Pstream::scatter(iMags_);
}


void Foam::functionObjects::STDMD::calcModes()
{
    Log<< tab << "# " << name() << ": Computing modes #" << endl;

    // Resize the number of output variables to nModes_, if need be
    if (nModes_ < iMags_.size())
    {
        iMags_.resize(nModes_);
    }

    // Compute and write out modes one by one
    label oModeI = 0;
    const RMatrix temp(QzUH_*RxInv_);

    for (const auto& i : iMags_)
    {
        List<complex> mode(temp.m(), Zero);

        // mode = temp*oEVecs_.subColumn(i);
        for (label p = 0; p < mode.size(); ++p)
        {
            for (label q = 0; q < oEVecs_.m(); ++q)
            {
                mode[p] += temp(p, q)*oEVecs_(q, i);
            }
        }

        volVectorField modeReal
        (
            IOobject
            (
                "modeReal" + std::to_string(oModeI) + fieldName_ + "_" + name(),
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector(dimless, Zero)
        );

        volVectorField modeImag
        (
            IOobject
            (
                "modeImag" + std::to_string(oModeI) + fieldName_ + "_" + name(),
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector(dimless, Zero)
        );

        label s = 0;
        for (label k = 0; k < nComps_; ++k)
        {
            for (label j = 0; j < modeReal.size(); ++j)
            {
                const complex& m = mode[s];
                modeReal[j].replace(k, m.real());
                modeImag[j].replace(k, m.imag());
                ++s;
            }
        }

        modeReal.write();
        modeImag.write();

        ++oModeI;
    }
}


Foam::scalar Foam::functionObjects::STDMD::sorter
(
    const List<scalar>& weight,
    const complex& amplitude,
    const complex& eval,
    const scalar modeNorm
) const
{
    // Omit eigenvalues with very large or very small magnitudes
    if (!(mag(eval) < GREAT && mag(eval) > VSMALL))
    {
        Info<< "  Returning zero magnitude for mag(eval) = " << mag(eval)
            << endl;

        return 0.0;
    }

    // Omit eigenvalue-STDMD step combinations that pose a risk of overflow
    if (mag(eval)*currIndex_ > sortLimiter_)
    {
        Info<< "  Returning zero magnitude for"
            << " mag(eval) = " << mag(eval)
            << " currIndex = " << currIndex_
            << " sortLimiter = " << sortLimiter_
            << endl;

        return 0.0;
    }

    scalar magnitude = 0;

    for (label j = 0; j < currIndex_; ++j)
    {
        magnitude += weight[j]*modeNorm*mag(amplitude*pow(eval, j + 1));
    }

    return magnitude;
}


void Foam::functionObjects::STDMD::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "STDMD output");
    writeCommented(os, "Frequency");
    writeTabbed(os, "Magnitude");
    writeTabbed(os, "Amplitude (real)");
    writeTabbed(os, "Amplitude (imag)");
    writeTabbed(os, "Eigenvalue (real)");
    writeTabbed(os, "Eigenvalue (imag)");
    os  << endl;
}


void Foam::functionObjects::STDMD::filterOutput()
{
    Log<< tab << "# " << name() << ": Filtering text output #" << endl;

    if (Pstream::master())
    {
        // Filter objects according to iMags
        filterIndexed(oEVals_, iMags_);
        filterIndexed(oEVecs_, iMags_);
        filterIndexed(oFreqs_, iMags_);
        filterIndexed(oAmps_, iMags_);
        filterIndexed(oMags_, iMags_);

        // Clip objects if need be (assuming objects have the same len)
        if (nModes_ < oFreqs_.size())
        {
            oEVals_.resize(nModes_);
            oEVecs_.resize(oEVecs_.m(), nModes_);
            oFreqs_.resize(nModes_);
            oAmps_.resize(nModes_);
            oMags_.resize(nModes_);
        }
    }
}


void Foam::functionObjects::STDMD::writeOutput(OFstream& os) const
{
    Log<< tab << "# " << name() << ": Writing text output #" << endl;

    writeFileHeader(os);

    for (const label i : labelRange(oFreqs_.size()))
    {
        os  << oFreqs_[i] << tab
            << oMags_[i] << tab
            << oAmps_[i].real() << tab
            << oAmps_[i].imag() << tab
            << oEVals_[i].real() << tab
            << oEVals_[i].imag() << endl;
    }
}


void Foam::functionObjects::STDMD::calcOutput()
{
    Log<< tab << "# " << name() << ":"
        << " Starts output processing for field = " << fieldName_ << " #"
        << endl;

    // Move upper and lower halves of Qz_ to new containers
    QzUH_ = Qz_.subMatrix(0, 0, nSnap_);
    QzLH_ = Qz_.subMatrix(nSnap_, 0, nSnap_);
    Qz_.clear();

    calcAp();

    QzLH_.clear();
    Gz_.clear();

    calcEigen();

    filterEVals();

    Ap_.clear();

    calcFreqs();

    calcFreqI();

    calcAmps();

    calcMags();

    calcMagI();

    calcModes();

    // Write out unfiltered data
    if (Pstream::master() && writeToFile_)
    {
        autoPtr<OFstream> osPtr =
            createFile("uSTDMD" + fieldName_, mesh_.time().timeOutputValue());
        OFstream& os = osPtr.ref();

        writeOutput(os);
    }

    filterOutput();

    // Write out filtered data
    if (Pstream::master() && writeToFile_)
    {
        autoPtr<OFstream> osPtr =
            createFile("STDMD" + fieldName_, mesh_.time().timeOutputValue());
        OFstream& os = osPtr.ref();

        writeOutput(os);
    }

    Log<< tab << "# " << name() << ":"
        << " Ends output processing for field = " << fieldName_ << " #"
        << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::STDMD::STDMD
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict, false),
    modeSorter_
    (
        modeSorterTypeNames.getOrDefault
        (
            "modeSorter",
            dict,
            modeSorterType::FIRST_SNAPSHOT
        )
    ),
    fieldName_(dict.get<word>("field")),
    initialised_(false),
    testEigen_(false),
    dumpEigen_(false),
    nModes_
    (
        dict.getCheckOrDefault<label>
        (
            "nModes",
            pTraits<label>::max,
            labelMinMax::ge(1)
        )
    ),
    maxRank_
    (
        dict.getCheckOrDefault<label>
        (
            "maxRank",
            pTraits<label>::max,
            labelMinMax::ge(1)
        )
    ),
    nGramSchmidt_
    (
        dict.getCheckOrDefault<label>
        (
            "nGramSchmidt",
            5,
            labelMinMax::ge(1)
        )
    ),
    fMin_
    (
        dict.getCheckOrDefault<label>
        (
            "fMin",
            0,
            labelMinMax::ge(0)
        )
    ),
    fMax_
    (
        dict.getCheckOrDefault<label>
        (
            "fMax",
            pTraits<label>::max,
            labelMinMax::ge(fMin_ + 1)
        )
    ),
    nComps_(Zero),
    nSnap_(Zero),
    currIndex_(Zero),
    sortLimiter_(500),
    absTol_(1e-2),
    relTol_(1e-3),
    minBasis_(Zero),
    dt_(Zero),
    minMagEVal_(Zero),
    zNorm_(Zero),
    ezNorm_(Zero),
    z_(),
    ez_(),
    X1_(),
    Qz_(),
    Gz_(),
    QzUH_(),
    QzLH_(),
    RxInv_(),
    Ap_(),
    oEVecs_(),
    oEVals_(),
    oAmps_(),
    oFreqs_(),
    iFreqs_(),
    oMags_(),
    iMags_()
{
    Info<< nl
        << "    ### This STDMD release is the beta release ###" << nl
        << "    Therefore small-to-medium changes in input/output interfaces "
        << "and internal structures should be expected in the next versions."
        << endl;

    // Check if varying or fixed time-step computation
    if (runTime.isAdjustTimeStep())
    {
        WarningInFunction
            << "  # STDMD: Viable only for fixed time-step computations. #"
            << endl;
    }

    // Check if mesh topology changing
    if (mesh_.topoChanging())
    {
        FatalErrorInFunction
            << "  # STDMD: Available only for non-changing mesh topology. #"
            << abort(FatalError);
    }

    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::STDMD::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    testEigen_ = dict.getOrDefault<bool>("testEigen", false);

    if (testEigen_)
    {
        dumpEigen_ = dict.getOrDefault<bool>("dumpEigen", false);
    }

    sortLimiter_ =
        dict.getCheckOrDefault<scalar>
        (
            "sortLimiter",
            500,
            scalarMinMax::ge(1)
        );

    absTol_ =
        dict.getCheckOrDefault<scalar>
        (
            "absTol",
            1e-2,
            scalarMinMax::ge(0)
        );

    relTol_ =
        dict.getCheckOrDefault<scalar>
        (
            "relTol",
            1e-3,
            scalarMinMax::ge(0)
        );

    minBasis_ =
        dict.getCheckOrDefault<scalar>
        (
            "minBasis",
            1e-8,
            scalarMinMax::ge(0)
        );

    dt_ =
        dict.getCheckOrDefault<scalar>
        (
            "stdmdInterval",
            (
                dict.getCheck<label>("executeInterval", labelMinMax::ge(1))
               *mesh_.time().deltaT().value()
            ),
            scalarMinMax::ge(0)
        );

    minMagEVal_ =
        dict.getCheckOrDefault<scalar>
        (
            "minEVal",
            1e-8,
            scalarMinMax::ge(0)
        );

    writeToFile_ = dict.getOrDefault<bool>("writeToFile", true);

    Info<< nl << name() << ":" << nl
        << tab << "# Input is read for field = "
        << fieldName_  << " #" << nl << endl;

    return true;
}


bool Foam::functionObjects::STDMD::execute()
{
    Log << type() << " " << name() << " execute:" << endl;

    // STDMD incremental orthonormal basis update (K:Fig. 15)
    snapshot();

    if (currIndex_ == 1)
    {
        initBasis();
    }

    if (1 < currIndex_)
    {
        GramSchmidt();

        // Check basis for z_ and, if necessary, expand Qz_ and Gz_
        zNorm_ = parnorm(z_);
        ezNorm_ = parnorm(ez_);

        if (minBasis_ < ezNorm_/zNorm_)
        {
            expandBasis();
        }

        updateGz();

        // Compress the orthonormal basis if required
        if (maxRank_ < Qz_.n())
        {
            compressBasis();
        }
    }
    ++currIndex_;

    Log<< tab << "# " << name() << ":"
        << " Execution index = " << currIndex_
        << " for field = " << fieldName_ << " #"
        << endl;

    return true;
}


bool Foam::functionObjects::STDMD::write()
{
    Log << type() << " " << name() << " write:" << endl;

    if (currIndex_ < 2)
    {
        WarningInFunction
            << "  # STDMD needs at least three snapshots to produce output #"
            << nl
            << "  # Only " << currIndex_ + 1 << " snapshots are available #"
            << nl
            << "  # Skipping STDMD output calculation and write #"
            << endl;

        return false;
    }

    // STDMD mode evaluation (K:Fig. 16)
    calcOutput();

    // Restart STDMD incremental orthonormal basis update
    initialised_ = false;

    mesh_.time().printExecutionTime(Info);

    return true;
}


// ************************************************************************* //
