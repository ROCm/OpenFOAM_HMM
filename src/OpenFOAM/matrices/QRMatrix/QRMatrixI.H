/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MatrixType>
inline void Foam::QRMatrix<MatrixType>::applyHouseholder
(
    MatrixType& A,
    const RMatrix& reflector,
    const label k
)
{
    applyLeftReflector(A, reflector, k, k);

    if (outputType_ == outputTypes::FULL_QR)
    {
        applyRightReflector(Q_, reflector, k);
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
template<class MatrixType>
inline Foam::RectangularMatrix<typename MatrixType::cmptType>
Foam::QRMatrix<MatrixType>::householderReflector
(
    RMatrix u
)
{
    #ifdef FULLDEBUG
    // Check if the given RectangularMatrix is effectively a column vector
    if (u.n() != 1)
    {
        FatalErrorInFunction
            << "Input matrix is not a column vector." << exit(FatalError);
    }
    #endif

    scalar magnitude(mag(u(0,0)));
    if (magnitude < VSMALL)
    {
        magnitude = SMALL;
        #if FULLDEBUG
        FatalErrorInFunction
            << "Almost zero leading elem in Householder matrix."
            << abort(FatalError);
        #endif
    }

    u(0,0) += u(0,0)/magnitude*u.columnNorm(0);

    scalar colNorm(u.columnNorm(0));
    if (colNorm < VSMALL)
    {
        colNorm = SMALL;
        #if FULLDEBUG
        FatalErrorInFunction
            << "Almost zero norm in the Householder matrix."
            << abort(FatalError);
        #endif
    }

    u /= cmptType(colNorm);

    return u;
}


template<class MatrixType>
inline void Foam::QRMatrix<MatrixType>::applyLeftReflector
(
    MatrixType& A,
    const RMatrix& reflector,
    const label k,
    const label k1
)
{
    //  const RMatrix& A0(A.subMatrix(k1, k));
    //  A0 -= (cmptType(2)*reflector)*(reflector & A0);

    for (label j = k; j < A.n(); ++j)
    {
        cmptType sum = Zero;
        for (label i = 0; i < reflector.m(); ++i)
        {
            sum += Detail::conj(reflector(i, 0))*A(i + k1, j);
        }

        sum *= cmptType(2);
        for (label i = 0; i < reflector.m(); ++i)
        {
            A(i + k1, j) -= reflector(i, 0)*sum;
        }
    }
}


template<class MatrixType>
inline void Foam::QRMatrix<MatrixType>::applyRightReflector
(
    MatrixType& A,
    const RMatrix& reflector,
    const label k
)
{
    // const RMatrix A0(A.subMatrix(0, k));
    // A0 -= ((A0*reflector)^(cmptType(2)*reflector))

    for (label i = 0; i < A.m(); ++i)
    {
        cmptType sum = Zero;
        for (label j = 0; j < reflector.m(); ++j)
        {
            sum += A(i, j + k)*reflector(j, 0);
        }

        sum *= cmptType(2);
        for (label j = 0; j < reflector.m(); ++j)
        {
            A(i, j + k) -= Detail::conj(reflector(j, 0))*sum;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MatrixType>
inline const MatrixType& Foam::QRMatrix<MatrixType>::Q() const
{
    const_cast<MatrixType&>(Q_).round();

    return Q_;
}


template<class MatrixType>
inline const MatrixType& Foam::QRMatrix<MatrixType>::R() const
{
    return R_;
}


template<class MatrixType>
inline const Foam::labelList& Foam::QRMatrix<MatrixType>::orderP() const
{
    return P_;
}


template<class MatrixType>
inline Foam::SquareMatrix<typename Foam::QRMatrix<MatrixType>::cmptType>
Foam::QRMatrix<MatrixType>::P() const
{
    SquareMatrix<cmptType> permMat(P_.size(), Zero);

    forAll(P_, jcol)
    {
        permMat(P_[jcol], jcol) = pTraits<cmptType>::one;
    }

    return permMat;
}


// ************************************************************************* //
