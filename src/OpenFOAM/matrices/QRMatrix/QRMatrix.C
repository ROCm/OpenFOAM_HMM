/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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

#include "QRMatrix.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MatrixType>
void Foam::QRMatrix<MatrixType>::qr
(
    MatrixType& A
)
{
    const label nIter = min(A.m() - 1, A.n());

    // Loop through all subcolumns of which the diagonal elem is the first elem
    for (label k = 0; k < nIter; ++k)
    {
        const RMatrix u(householderReflector(A.subColumn(k, k)));

        applyHouseholder(A, u, k);
    }

    if (outputType_ == outputTypes::REDUCED_R)
    {
        A.resize(A.n(), A.n());
    }
}


template<class MatrixType>
void Foam::QRMatrix<MatrixType>::qrPivot
(
    MatrixType& A
)
{
    const label nCols = A.n();
    const label nIter = min(A.m() - 1, A.n());

    // Initialise permutation vector, and column norm vector
    P_ = identity(nCols);

    // Initialise vector norms of each column of A
    List<scalar> colNorms(nCols);
    for (label k = 0; k < nCols; ++k)
    {
        colNorms[k] = A.columnNorm(k, true);
    }

    // Loop through all subcolumns of which the diagonal elem is the first elem
    for (label k = 0; k < nIter; ++k)
    {
        const labelRange colRange(k, nCols);
        const SubList<scalar> subColNorms(colNorms, colRange);

        // Column pivoting
        const label maxColNormi =
            std::distance
            (
                subColNorms.cbegin(),
                std::max_element(subColNorms.cbegin(), subColNorms.cend())
            );

        // Swap R_, P_ and colNorms_ according to pivot column if the current
        // column is not the max norm column by avoiding any column swap where
        // the leading elem is very small
        if (maxColNormi != 0 && SMALL < mag(A(k, k + maxColNormi)))
        {
            const RMatrix R1(A.subColumn(k));
            const RMatrix R2(A.subColumn(maxColNormi + k));
            A.subColumn(k) = R2;
            A.subColumn(maxColNormi + k) = R1;

            Swap(P_[k], P_[maxColNormi + k]);
            Swap(colNorms[k], colNorms[maxColNormi + k]);
        }

        {
            const RMatrix u(householderReflector(A.subColumn(k, k)));

            applyHouseholder(A, u, k);
        }

        // Update norms
        if (k < nIter - 1)
        {
            label q = k + 1;
            for (const auto& val : RMatrix(A.subRow(k, q)))
            {
                colNorms[q] -= magSqr(val);
                ++q;
            }
        }
    }

    if (outputType_ == outputTypes::REDUCED_R)
    {
        A.resize(A.n(), A.n());
    }
}


template<class MatrixType>
template<template<typename> class ListContainer>
void Foam::QRMatrix<MatrixType>::solvex
(
    ListContainer<cmptType>& x
) const
{
    const label n = R_.n();

    for (label i = n - 1; 0 <= i; --i)
    {
        cmptType sum = x[i];

        for (label j = i + 1; j < n; ++j)
        {
            sum -= x[j]*R_(i, j);
        }

        if (mag(R_(i, i)) < SMALL)
        {
            FatalErrorInFunction
                << "Back-substitution failed due to small diagonal"
                << abort(FatalError);
        }

        x[i] = sum/R_(i, i);
    }
}


template<class MatrixType>
template<template<typename> class ListContainer>
void Foam::QRMatrix<MatrixType>::solveImpl
(
    List<cmptType>& x,
    const ListContainer<cmptType>& source
) const
{
    // Assert (&x != &source) ?
    const label m = Q_.m();

    // x = Q_.T()*source;  i.e., Q_.Tmul(source)
    for (label i = 0; i < m; ++i)
    {
        x[i] = 0;
        for (label j = 0; j < m; ++j)
        {
            x[i] += Q_(j, i)*source[j];
        }
    }

    solvex(x);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MatrixType>
Foam::QRMatrix<MatrixType>::QRMatrix()
:
    outputType_(),
    storeMethod_(),
    colPivot_()
{}


template<class MatrixType>
Foam::QRMatrix<MatrixType>::QRMatrix
(
    const outputTypes outputType,
    const storeMethods storeMethod,
    const colPivoting colPivot
)
:
    outputType_(outputType),
    storeMethod_(storeMethod),
    colPivot_(colPivot),
    Q_(),
    R_(),
    P_()
{}


template<class MatrixType>
Foam::QRMatrix<MatrixType>::QRMatrix
(
    MatrixType& A,
    const outputTypes outputType,
    const storeMethods storeMethod,
    const colPivoting colPivot
)
:
    outputType_(outputType),
    storeMethod_(storeMethod),
    colPivot_(colPivot),
    Q_(),
    R_(),
    P_()
{
    decompose(A);
}


template<class MatrixType>
Foam::QRMatrix<MatrixType>::QRMatrix
(
    const MatrixType& A,
    const outputTypes outputType,
    const colPivoting colPivot
)
:
    outputType_(outputType),
    storeMethod_(storeMethods::OUT_OF_PLACE),
    colPivot_(colPivot),
    Q_(),
    R_(),
    P_()
{
    decompose(A);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class MatrixType>
void Foam::QRMatrix<MatrixType>::decompose
(
    MatrixType& A
)
{
    // Check whether settings and input are consistent for reduced QRMatrix
    if (A.m() <= A.n() && outputType_ == outputTypes::REDUCED_R)
    {
        outputType_ = outputTypes::FULL_R;

        #if FULLDEBUG
        WarningInFunction
            << "Reduced QR decomposition is by definition limited to matrices"
            << " wherein rows > columns, thus computing FULL decomposition."
            << nl;
        #endif
    }

    // Allocate resources for Q_, if need be
    if (outputType_ == outputTypes::FULL_QR)
    {
        Q_ = MatrixType({A.m(), A.m()}, I);
    }

    // Allocate resources for R_ and execute decomposition
    switch (storeMethod_)
    {
        case storeMethods::IN_PLACE:
        {
            if (colPivot_)
            {
                qrPivot(A);
            }
            else
            {
                qr(A);
            }
            break;
        }

        case storeMethods::OUT_OF_PLACE:
        {
            R_ = A;

            if (colPivot_)
            {
                qrPivot(R_);
            }
            else
            {
                qr(R_);
            }
            break;
        }
    }
}


template<class MatrixType>
void Foam::QRMatrix<MatrixType>::decompose
(
    const MatrixType& A
)
{
    if (storeMethod_ == storeMethods::IN_PLACE)
    {
        WarningInFunction
            << "const type qualifier invalidates storeMethods::IN_PLACE." << nl;
    }

    // Check whether settings and input are consistent for reduced QRMatrix
    if (A.m() <= A.n() && outputType_ == outputTypes::REDUCED_R)
    {
        outputType_ = outputTypes::FULL_R;

        #if FULLDEBUG
        WarningInFunction
            << "Reduced QR decomposition is by definition limited to matrices"
            << " wherein rows > columns, thus computing FULL decomposition."
            << nl;
        #endif
    }

    // Allocate resources for Q_, if need be
    if (outputType_ == outputTypes::FULL_QR)
    {
        Q_ = MatrixType({A.m(), A.m()}, I);
    }

    // Allocate resources for R_ and execute decomposition
    R_ = A;

    if (colPivot_)
    {
        qrPivot(R_);
    }
    else
    {
        qr(R_);
    }
}


template<class MatrixType>
void Foam::QRMatrix<MatrixType>::solve
(
    List<cmptType>& x,
    const UList<cmptType>& source
) const
{
    solveImpl(x, source);
}


template<class MatrixType>
template<class Addr>
void Foam::QRMatrix<MatrixType>::solve
(
    List<cmptType>& x,
    const IndirectListBase<cmptType, Addr>& source
) const
{
    solveImpl(x, source);
}


template<class MatrixType>
Foam::tmp<Foam::Field<typename MatrixType::cmptType>>
Foam::QRMatrix<MatrixType>::solve
(
    const UList<cmptType>& source
) const
{
    auto tresult(Q_.Tmul(source));

    solvex(tresult.ref());

    return tresult;
}


template<class MatrixType>
template<class Addr>
Foam::tmp<Foam::Field<typename MatrixType::cmptType>>
Foam::QRMatrix<MatrixType>::solve
(
    const IndirectListBase<cmptType, Addr>& source
) const
{
    auto tresult(Q_.Tmul(source));

    solvex(tresult.ref());

    return tresult;
}


template<class MatrixType>
typename Foam::QRMatrix<MatrixType>::SMatrix
Foam::QRMatrix<MatrixType>::inv() const
{
    const label m = Q_.m();

    Field<cmptType> x(m);
    SMatrix inv(m);

    for (label i = 0; i < m; ++i)
    {
        for (label j = 0; j < m; ++j)
        {
            x[j] = Q_(i, j);
        }
        solvex(x);
        inv.subColumn(i) = x;
    }

    return inv;
}


template<class MatrixType>
Foam::RectangularMatrix<typename MatrixType::cmptType>
Foam::QRMatrix<MatrixType>::backSubstitution
(
    const SMatrix& A,
    const RMatrix& rhs
)
{
    const label mRows = A.m();
    const label nCols = A.n();
    const label pCols = rhs.n();

    #ifdef FULLDEBUG
    {
        const label qRows = rhs.m();
        if (mRows != qRows)
        {
            FatalErrorInFunction
                << "Linear system is not solvable since the number of rows of"
                << "A and rhs are not equal:" << tab << mRows << "vs." << qRows
                << abort(FatalError);
        }

        const List<cmptType> diag(A.diag());
        for (const auto& val : diag)
        {
            if (mag(val) < SMALL)
            {
                WarningInFunction
                    << "SMALL-valued diagonal elem in back-substitution." << nl;
            }
        }
    }
    #endif

    RMatrix X({nCols, pCols}, Zero);

    for (label i = 0; i < pCols; ++i)
    {
        for (label j = mRows - 1; -1 < j; --j)
        {
            cmptType alpha(rhs(j, i));

            for (label k = j + 1; k < mRows; ++k)
            {
                alpha -= X(k, i)*A(j, k);
            }

            X(j, i) = alpha/A(j, j);
        }
    }

    return X;
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class MatrixType>
MatrixType Foam::pinv
(
    const MatrixType& A,
    const scalar tolerance
)
{
    scalar tol = tolerance;
    typedef typename MatrixType::cmptType cmptType;

    if (A.empty())
    {
        FatalErrorInFunction
            << "Empty matrix found."
            << abort(FatalError);
    }

    if (A.size() == 1)
    {
        if (A(0,0) == cmptType(0))
        {
            return MatrixType({1, 1}, cmptType(0));
        }

        return MatrixType({1, 1}, cmptType(1)/(A(0,0) + cmptType(VSMALL)));
    }

    QRMatrix<MatrixType> QRM
    (
        A,
        QRMatrix<MatrixType>::outputTypes::FULL_QR,
        QRMatrix<MatrixType>::colPivoting::TRUE
    );
    const MatrixType& R(QRM.R());
    const MatrixType& Q(QRM.Q());

    // R1 (KP:p. 648)
    // Find the first diagonal element index with (almost) zero value in R
    label firstZeroElemi = 0;
    label i = 0;
    while (i < 2)
    {
        const List<cmptType> diag(R.diag());

        auto lessThan = [=](const cmptType& x) { return tol > mag(x); };

        firstZeroElemi =
            std::distance
            (
                diag.cbegin(),
                std::find_if(diag.cbegin(), diag.cend(), lessThan)
            );

        if (firstZeroElemi == 0)
        {
            if (i == 0)
            {
                WarningInFunction
                    << "The largest diagonal element magnitude is nearly zero. "
                    << "Tightening the tolerance."
                    << endl;

                tol = 1e-13;
                ++i;
            }
            else
            {
                WarningInFunction
                    << "The largest diagonal element magnitude is nearly zero. "
                    << "Returning a zero matrix."
                    << endl;
                ++i;

                return MatrixType({A.n(), A.m()}, Zero);
            }
        }
        else
        {
            i += 2;
        }
    }

    // Exclude the first (almost) zero diagonal row and the rows below
    // since R diagonal is already descending due to the QR column pivoting
    const RectangularMatrix<cmptType> R1(R.subMatrix(0, 0, firstZeroElemi));

    // R2 (KP:p. 648)
    if (R1.n() < R1.m())
    {
        const SquareMatrix<cmptType> C(R1&R1);

        QRMatrix<SquareMatrix<cmptType>> QRSolve
        (
            C,
            QRMatrix<SquareMatrix<cmptType>>::outputTypes::FULL_QR
        );

        RectangularMatrix<cmptType> R2
        (
            QRSolve.backSubstitution
            (
                QRSolve.R(),
                QRSolve.Q() & (R1.T())
            )
        );

        // R3 (KP:p. 648)
        R2.resize(R.m(), R.n());

        return MatrixType((QRM.P()^R2)^Q);
    }
    else
    {
        const SquareMatrix<cmptType> C(R1^R1);

        QRMatrix<SquareMatrix<cmptType>> QRSolve
        (
            C,
            QRMatrix<SquareMatrix<cmptType>>::outputTypes::FULL_QR
        );

        RectangularMatrix<cmptType> R2
        (
            QRSolve.backSubstitution
            (
                QRSolve.R(),
                QRSolve.Q() & R1
            )
        );

        // R3
        R2.resize(R.m(), R.n());

        return MatrixType((QRM.P()^R2)^Q);
    }
}


// ************************************************************************* //
