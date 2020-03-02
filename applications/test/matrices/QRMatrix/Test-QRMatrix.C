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

#include "MatrixTools.H"
#include "QRMatrix.H"
#include "Random.H"
#include "IOmanip.H"

using namespace Foam;
using namespace Foam::MatrixTools;

#define equal MatrixTools::equal
#define RUNALL true
const bool verbose = true;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void horizontalLine()
{
    Info<< "+---------+---------+---------+---------+---------+" << nl;
}


// Create random scalar-type matrix
template<class MatrixType>
typename std::enable_if
<
    !std::is_same<complex, typename MatrixType::cmptType>::value,
    MatrixType
>::type makeRandomMatrix
(
    const labelPair& dims,
    Random& rndGen
)
{
    MatrixType mat(dims);

    std::generate
    (
        mat.begin(),
        mat.end(),
        [&]{ return rndGen.GaussNormal<typename MatrixType::cmptType>(); }
    );

    return mat;
}


// Create random complex-type matrix
template<class MatrixType>
typename std::enable_if
<
    std::is_same<complex, typename MatrixType::cmptType>::value,
    MatrixType
>::type makeRandomMatrix
(
    const labelPair& dims,
    Random& rndGen
)
{
    MatrixType mat(dims);

    std::generate
    (
        mat.begin(),
        mat.end(),
        [&]
        {
            return complex
            (
                rndGen.GaussNormal<scalar>(),
                rndGen.GaussNormal<scalar>()
            );
        }
    );

    return mat;
}


// Print OpenFOAM matrix in NumPy format
template<class MatrixType>
void InfoNumPyFormat
(
    const MatrixType& mat,
    const word objName
)
{
    Info<< objName << ": " << mat.m() << "x" << mat.n() << nl;
    for (label m = 0; m < mat.m(); ++m)
    {
        Info<< "[";
        for (label n = 0; n < mat.n(); ++n)
        {
            if (n == mat.n() - 1)
            {
                Info<< mat(m,n);
            }
            else
            {
                Info<< mat(m,n) << ",";
            }
        }
        if (m == mat.m() - 1)
        {
            Info << "]" << nl;
        }
        else
        {
            Info << "]," << nl;
        }
    }
}


// Returns true if two scalars are equal within a given tolerance, and v.v.
bool isEqual
(
    const scalar a,
    const scalar b,
    const bool verbose = false,
    const scalar relTol = 1e-5,
    const scalar absTol = 1e-8
)
{
    if ((absTol + relTol*mag(b)) < mag(a - b))
    {
        if (verbose)
        {
            Info<< "Elements are not close in terms of tolerances:"
                << nl << a << tab << b << nl;
        }
        return false;
    }

    if (verbose)
    {
        Info<< "All elems are the same within the tolerances" << nl;
    }

    return true;
}


// Prints true if a given matrix is empty, and v.v.
template<class MatrixType>
void isEmpty
(
    const MatrixType& A,
    const word objName
)
{
    Info<< "Empty " << objName << " = ";
    if (A.empty())
    {
        Info<< "true" << nl;
    }
    else
    {
        Info<< "false" << nl;
    }
}


// Checks if Q matrix is a unitary matrix
template<class MatrixType>
void cross_check_QRMatrix
(
    const MatrixType& Aorig,
    const MatrixType& Q
)
{
    InfoNumPyFormat(Q, "Q");

    // mathworld.wolfram.com/UnitaryMatrix.html (Retrieved:18-06-19)
    Info<< nl << "# Q.T()*Q = Q*Q.T() = I:" << nl;
    const MatrixType QTQ(Q&Q);
    const MatrixType QQT(Q^Q);
    MatrixType IMatrix({Q.m(), Q.n()}, Zero);
    for (label i = 0; i < IMatrix.n(); ++i)
    {
        IMatrix(i,i) = pTraits<typename MatrixType::cmptType>::one;
    }
    equal(QTQ, IMatrix, verbose);
    equal(QQT, IMatrix, verbose);
    equal(QTQ, QQT, verbose);
}


// Checks if R matrix is an upper-triangular matrix
template<class MatrixType>
void cross_check_QRMatrix
(
    const MatrixType& R
)
{
    InfoNumPyFormat(R, "R");

    // mathworld.wolfram.com/UpperTriangularMatrix.html (Retrieved:18-06-19)
    Info<< nl << "# R(i, j) = 0 for i > j:" << nl;
    for (label i = 0; i < R.m(); ++i)
    {
        for (label j = 0; j < R.n(); ++j)
        {
            if (j < i)
            {
                isEqual(mag(R(i, j)), 0.0, verbose);
            }
        }
    }
}


// Checks if the given matrix can be reconstructed by Q and R matrices
template<class MatrixType>
void cross_check_QRMatrix
(
    const MatrixType& Aorig,
    const MatrixType& Q,
    const MatrixType& R
)
{
    InfoNumPyFormat(Aorig, "Aorig");

    cross_check_QRMatrix(Aorig, Q);

    cross_check_QRMatrix(R);

    // mathworld.wolfram.com/QRDecomposition.html (Retrieved:18-06-19)
    Info<< nl << "# Q*R = A:" << nl;
    const MatrixType AReconstruct(Q*R);
    equal(Aorig, AReconstruct, verbose);
}


// Checks if R matrix is an upper-triangular matrix, and obeys column pivoting
template<class MatrixType>
void cross_check_QRMatrix
(
    const MatrixType& Aorig,
    const labelList& orderP,
    const SquareMatrix<typename MatrixType::cmptType>& P,
    const MatrixType& R
)
{
    InfoNumPyFormat(Aorig, "Aorig");

    InfoNumPyFormat(P, "P");

    cross_check_QRMatrix(R);

    // w.wiki/54m (Retrieved:20-06-19)
    Info<< nl << "# Diag elems of R non-increasing:"
        << "|R11| >= |R22| >= ... >= |Rnn|:" << nl;
    const List<scalar> diag0(mag(R.diag()));
    List<scalar> diag1(diag0);
    Foam::sort(diag1, std::greater<scalar>());
    forAll(diag0, i)
    {
        isEqual(diag0[i], diag1[i], verbose);
    }
}


// Checks if the given matrix can be reconstructed by column-pivoted Q and R
template<class MatrixType>
void cross_check_QRMatrix
(
    const MatrixType& Aorig,
    const MatrixType& Q,
    const labelList& orderP,
    const SquareMatrix<typename MatrixType::cmptType>& P,
    const MatrixType& R
)
{
    cross_check_QRMatrix(Aorig, orderP, P, R);
    cross_check_QRMatrix(Aorig, Q);

    // w.wiki/54m (Retrieved:20-06-19)
    Info<< nl << "# Q*R = A*P:" << nl;
    const MatrixType AReconstruct(Q*R);
    const MatrixType AP(Aorig*P);
    equal(AP, AReconstruct, verbose);
}


// Checks each constructor of QRMatrix type
template<class MatrixType>
void verification_QRMatrix
(
    MatrixType& A
)
{
    typedef SquareMatrix<typename MatrixType::cmptType> SMatrix;

    // Create working copies of matrix A
    const MatrixType Aorig(A);

    // QRMatrix Constructors
    #if (0 | RUNALL)
    {
        QRMatrix<MatrixType> QRNull;
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# FULL_R | OUT_OF_PLACE" << nl;
        QRMatrix<MatrixType> QRM
        (
            QRMatrix<MatrixType>::outputTypes::FULL_R,
            QRMatrix<MatrixType>::storeMethods::OUT_OF_PLACE
        );
        QRM.decompose(A);
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        isEmpty(Q, "Q");
        cross_check_QRMatrix(R);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# FULL_R | IN_PLACE" << nl;
        QRMatrix<MatrixType> QRM
        (
            QRMatrix<MatrixType>::outputTypes::FULL_R,
            QRMatrix<MatrixType>::storeMethods::IN_PLACE
        );
        MatrixType A0(A);
        QRM.decompose(A0);
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        isEmpty(Q, "Q");
        isEmpty(R, "R");
        cross_check_QRMatrix(A0);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# FULL_QR | OUT_OF_PLACE" << nl;
        QRMatrix<MatrixType> QRM
        (
            QRMatrix<MatrixType>::outputTypes::FULL_QR,
            QRMatrix<MatrixType>::storeMethods::OUT_OF_PLACE
        );
        QRM.decompose(A);
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        cross_check_QRMatrix(Aorig, Q, R);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# FULL_QR | IN_PLACE" << nl;
        QRMatrix<MatrixType> QRM
        (
            QRMatrix<MatrixType>::outputTypes::FULL_QR,
            QRMatrix<MatrixType>::storeMethods::IN_PLACE
        );
        MatrixType A0(A);
        QRM.decompose(A0);
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        isEmpty(R, "R");
        cross_check_QRMatrix(Aorig, Q, A0);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# FULL_R | OUT_OF_PLACE | colPivot = on" << nl;
        QRMatrix<MatrixType> QRM
        (
            QRMatrix<MatrixType>::outputTypes::FULL_R,
            QRMatrix<MatrixType>::storeMethods::OUT_OF_PLACE,
            QRMatrix<MatrixType>::colPivoting::TRUE
        );
        QRM.decompose(A);
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        const labelList& PList(QRM.orderP());
        const SMatrix P(QRM.P());
        isEmpty(Q, "Q");
        cross_check_QRMatrix(Aorig, PList, P, R);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# FULL_R | IN_PLACE | colPivot = on" << nl;
        QRMatrix<MatrixType> QRM
        (
            QRMatrix<MatrixType>::outputTypes::FULL_R,
            QRMatrix<MatrixType>::storeMethods::IN_PLACE,
            QRMatrix<MatrixType>::colPivoting::TRUE
        );
        MatrixType A0(A);
        QRM.decompose(A0);
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        const labelList& PList(QRM.orderP());
        const SMatrix P(QRM.P());
        isEmpty(Q, "Q");
        isEmpty(R, "R");
        cross_check_QRMatrix(Aorig, PList, P, A0);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# FULL_QR | OUT_OF_PLACE | colPivot = on" << nl;
        QRMatrix<MatrixType> QRM
        (
            QRMatrix<MatrixType>::outputTypes::FULL_QR,
            QRMatrix<MatrixType>::storeMethods::OUT_OF_PLACE,
            QRMatrix<MatrixType>::colPivoting::TRUE
        );
        QRM.decompose(A);
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        const labelList& PList(QRM.orderP());
        const SMatrix P(QRM.P());
        cross_check_QRMatrix(Aorig, Q, PList, P, R);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# FULL_QR | IN_PLACE | colPivot = on" << nl;
        QRMatrix<MatrixType> QRM
        (
            QRMatrix<MatrixType>::outputTypes::FULL_QR,
            QRMatrix<MatrixType>::storeMethods::IN_PLACE,
            QRMatrix<MatrixType>::colPivoting::TRUE
        );
        MatrixType A0(A);
        QRM.decompose(A0);
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        const labelList& PList(QRM.orderP());
        const SMatrix P(QRM.P());
        isEmpty(R, "R");
        cross_check_QRMatrix(Aorig, Q, PList, P, A0);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# A | FULL_R | OUT_OF_PLACE | colPivot = off" << nl;
        MatrixType A0(A);
        QRMatrix<MatrixType> QRM
        (
            A0,
            QRMatrix<MatrixType>::outputTypes::FULL_R,
            QRMatrix<MatrixType>::storeMethods::OUT_OF_PLACE,
            QRMatrix<MatrixType>::colPivoting::FALSE
        );
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        isEmpty(Q, "Q");
        cross_check_QRMatrix(R);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# A | FULL_R | IN_PLACE | colPivot = off" << nl;
        MatrixType A0(A);
        QRMatrix<MatrixType> QRM
        (
            A0,
            QRMatrix<MatrixType>::outputTypes::FULL_R,
            QRMatrix<MatrixType>::storeMethods::IN_PLACE,
            QRMatrix<MatrixType>::colPivoting::FALSE
        );
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        isEmpty(Q, "Q");
        isEmpty(R, "R");
        cross_check_QRMatrix(A0);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# A | FULL_QR | OUT_OF_PLACE | colPivot = off" << nl;
        MatrixType A0(A);
        QRMatrix<MatrixType> QRM
        (
            A0,
            QRMatrix<MatrixType>::outputTypes::FULL_QR,
            QRMatrix<MatrixType>::storeMethods::OUT_OF_PLACE,
            QRMatrix<MatrixType>::colPivoting::FALSE
        );
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        cross_check_QRMatrix(Aorig, Q, R);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# A | FULL_QR | IN_PLACE | colPivot = off" << nl;
        MatrixType A0(A);
        QRMatrix<MatrixType> QRM
        (
            A0,
            QRMatrix<MatrixType>::outputTypes::FULL_QR,
            QRMatrix<MatrixType>::storeMethods::IN_PLACE,
            QRMatrix<MatrixType>::colPivoting::FALSE
        );
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        isEmpty(R, "R");
        cross_check_QRMatrix(Aorig, Q, A0);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# A | FULL_R | OUT_OF_PLACE | colPivot = on" << nl;
        MatrixType A0(A);
        QRMatrix<MatrixType> QRM
        (
            A0,
            QRMatrix<MatrixType>::outputTypes::FULL_R,
            QRMatrix<MatrixType>::storeMethods::OUT_OF_PLACE,
            QRMatrix<MatrixType>::colPivoting::TRUE
        );
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        const labelList& PList(QRM.orderP());
        const SMatrix P(QRM.P());
        isEmpty(Q, "Q");
        cross_check_QRMatrix(Aorig, PList, P, R);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# A | FULL_R | IN_PLACE | colPivot = on" << nl;
        MatrixType A0(A);
        QRMatrix<MatrixType> QRM
        (
            A0,
            QRMatrix<MatrixType>::outputTypes::FULL_R,
            QRMatrix<MatrixType>::storeMethods::IN_PLACE,
            QRMatrix<MatrixType>::colPivoting::TRUE
        );
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        const labelList& PList(QRM.orderP());
        const SMatrix P(QRM.P());
        isEmpty(Q, "Q");
        isEmpty(R, "R");
        cross_check_QRMatrix(Aorig, PList, P, A0);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# A | FULL_QR | OUT_OF_PLACE | colPivot = on" << nl;
        MatrixType A0(A);
        QRMatrix<MatrixType> QRM
        (
            A0,
            QRMatrix<MatrixType>::outputTypes::FULL_QR,
            QRMatrix<MatrixType>::storeMethods::OUT_OF_PLACE,
            QRMatrix<MatrixType>::colPivoting::TRUE
        );
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        const labelList& PList(QRM.orderP());
        const SMatrix P(QRM.P());
        cross_check_QRMatrix(Aorig, Q, PList, P, R);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# A | FULL_QR | IN_PLACE | colPivot = on" << nl;
        MatrixType A0(A);
        QRMatrix<MatrixType> QRM
        (
            A0,
            QRMatrix<MatrixType>::outputTypes::FULL_QR,
            QRMatrix<MatrixType>::storeMethods::IN_PLACE,
            QRMatrix<MatrixType>::colPivoting::TRUE
        );
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        const labelList& PList(QRM.orderP());
        const SMatrix P(QRM.P());
        isEmpty(R, "R");
        cross_check_QRMatrix(Aorig, Q, PList, P, A0);
    }
    #endif

    // constructors with the const type specifier
    #if (0 | RUNALL)
    {
        Info<< "# const A | FULL_R | colPivot = off" << nl;
        const MatrixType A0(A);
        QRMatrix<MatrixType> QRM
        (
            A0,
            QRMatrix<MatrixType>::outputTypes::FULL_R,
            QRMatrix<MatrixType>::colPivoting::FALSE
        );
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        isEmpty(Q, "Q");
        cross_check_QRMatrix(R);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# const A | FULL_QR | colPivot = off" << nl;
        const MatrixType A0(A);
        QRMatrix<MatrixType> QRM
        (
            A0,
            QRMatrix<MatrixType>::outputTypes::FULL_QR,
            QRMatrix<MatrixType>::colPivoting::FALSE
        );
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        cross_check_QRMatrix(Aorig, Q, R);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# const A | FULL_R | colPivot = on" << nl;
        const MatrixType A0(A);
        QRMatrix<MatrixType> QRM
        (
            A0,
            QRMatrix<MatrixType>::outputTypes::FULL_R,
            QRMatrix<MatrixType>::colPivoting::TRUE
        );
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        const labelList& PList(QRM.orderP());
        const SMatrix P(QRM.P());
        isEmpty(Q, "Q");
        cross_check_QRMatrix(Aorig, PList, P, R);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# const A | FULL_QR | colPivot = on" << nl;
        const MatrixType A0(A);
        QRMatrix<MatrixType> QRM
        (
            A0,
            QRMatrix<MatrixType>::outputTypes::FULL_QR,
            QRMatrix<MatrixType>::colPivoting::TRUE
        );
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        const labelList& PList(QRM.orderP());
        const SMatrix P(QRM.P());
        cross_check_QRMatrix(Aorig, Q, PList, P, R);
    }
    #endif

    //
    #if (0 | RUNALL)
    {
        Info<< "# REDUCED_R | OUT_OF_PLACE" << nl;
        QRMatrix<MatrixType> QRM
        (
            QRMatrix<MatrixType>::outputTypes::REDUCED_R,
            QRMatrix<MatrixType>::storeMethods::OUT_OF_PLACE
        );
        QRM.decompose(A);
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        isEmpty(Q, "Q");
        cross_check_QRMatrix(R);

        // check if full R and reduced R give the same answer
        if (A.n() < A.m())
        {
            QRMatrix<MatrixType> fullQRM
            (
                QRMatrix<MatrixType>::outputTypes::FULL_QR,
                QRMatrix<MatrixType>::storeMethods::OUT_OF_PLACE
            );
            fullQRM.decompose(A);
            const MatrixType& fullR(fullQRM.R());
            MatrixType fullR2(fullR.subMatrix(0,0,R.m()));
            equal(fullR2, R, verbose);
        }
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# REDUCED_R | IN_PLACE" << nl;
        QRMatrix<MatrixType> QRM
        (
            QRMatrix<MatrixType>::outputTypes::REDUCED_R,
            QRMatrix<MatrixType>::storeMethods::IN_PLACE
        );
        MatrixType A0(A);
        QRM.decompose(A0);
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        isEmpty(Q, "Q");
        isEmpty(R, "R");
        cross_check_QRMatrix(A0);

        // check if full R and reduced R give the same answer
        if (A.n() < A.m())
        {
            QRMatrix<MatrixType> fullQRM
            (
                QRMatrix<MatrixType>::outputTypes::FULL_QR,
                QRMatrix<MatrixType>::storeMethods::IN_PLACE
            );
            MatrixType A1(A);
            fullQRM.decompose(A1);
            MatrixType fullR(A1.subMatrix(0,0,A0.m()));
            equal(fullR, A0, verbose);
        }
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# REDUCED_R | OUT_OF_PLACE | colPivot = on" << nl;
        QRMatrix<MatrixType> QRM
        (
            QRMatrix<MatrixType>::outputTypes::REDUCED_R,
            QRMatrix<MatrixType>::storeMethods::OUT_OF_PLACE,
            QRMatrix<MatrixType>::colPivoting::TRUE
        );
        QRM.decompose(A);
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        const labelList& PList(QRM.orderP());
        const SMatrix P(QRM.P());
        isEmpty(Q, "Q");
        cross_check_QRMatrix(Aorig, PList, P, R);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# REDUCED_R | IN_PLACE | colPivot = on" << nl;
        QRMatrix<MatrixType> QRM
        (
            QRMatrix<MatrixType>::outputTypes::REDUCED_R,
            QRMatrix<MatrixType>::storeMethods::IN_PLACE,
            QRMatrix<MatrixType>::colPivoting::TRUE
        );
        MatrixType A0(A);
        QRM.decompose(A0);
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        const labelList& PList(QRM.orderP());
        const SMatrix P(QRM.P());
        isEmpty(Q, "Q");
        isEmpty(R, "R");
        cross_check_QRMatrix(Aorig, PList, P, A0);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# A | REDUCED_R | OUT_OF_PLACE | colPivot = off" << nl;
        MatrixType A0(A);
        QRMatrix<MatrixType> QRM
        (
            A0,
            QRMatrix<MatrixType>::outputTypes::REDUCED_R,
            QRMatrix<MatrixType>::storeMethods::OUT_OF_PLACE,
            QRMatrix<MatrixType>::colPivoting::FALSE
        );
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        isEmpty(Q, "Q");
        cross_check_QRMatrix(R);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# A | REDUCED_R | IN_PLACE | colPivot = off" << nl;
        MatrixType A0(A);
        QRMatrix<MatrixType> QRM
        (
            A0,
            QRMatrix<MatrixType>::outputTypes::REDUCED_R,
            QRMatrix<MatrixType>::storeMethods::IN_PLACE,
            QRMatrix<MatrixType>::colPivoting::FALSE
        );
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        isEmpty(Q, "Q");
        isEmpty(R, "R");
        cross_check_QRMatrix(A0);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# A | REDUCED_R | OUT_OF_PLACE | colPivot = on" << nl;
        MatrixType A0(A);
        QRMatrix<MatrixType> QRM
        (
            A0,
            QRMatrix<MatrixType>::outputTypes::REDUCED_R,
            QRMatrix<MatrixType>::storeMethods::OUT_OF_PLACE,
            QRMatrix<MatrixType>::colPivoting::TRUE
        );
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        const labelList& PList(QRM.orderP());
        const SMatrix P(QRM.P());
        isEmpty(Q, "Q");
        cross_check_QRMatrix(Aorig, PList, P, R);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# A | REDUCED_R | IN_PLACE | colPivot = on" << nl;
        MatrixType A0(A);
        QRMatrix<MatrixType> QRM
        (
            A0,
            QRMatrix<MatrixType>::outputTypes::REDUCED_R,
            QRMatrix<MatrixType>::storeMethods::IN_PLACE,
            QRMatrix<MatrixType>::colPivoting::TRUE
        );
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        const labelList& PList(QRM.orderP());
        const SMatrix P(QRM.P());
        isEmpty(Q, "Q");
        isEmpty(R, "R");
        cross_check_QRMatrix(Aorig, PList, P, A0);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# const A | REDUCED_R | colPivot = off" << nl;
        const MatrixType A0(A);
        QRMatrix<MatrixType> QRM
        (
            A0,
            QRMatrix<MatrixType>::outputTypes::REDUCED_R,
            QRMatrix<MatrixType>::colPivoting::FALSE
        );
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        isEmpty(Q, "Q");
        cross_check_QRMatrix(R);
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< "# const A | REDUCED_R | colPivot = on" << nl;
        const MatrixType A0(A);
        QRMatrix<MatrixType> QRM
        (
            A0,
            QRMatrix<MatrixType>::outputTypes::REDUCED_R,
            QRMatrix<MatrixType>::colPivoting::TRUE
        );
        const MatrixType& Q(QRM.Q());
        const MatrixType& R(QRM.R());
        const labelList& PList(QRM.orderP());
        const SMatrix P(QRM.P());
        isEmpty(Q, "Q");
        cross_check_QRMatrix(Aorig, PList, P, R);
    }
    #endif
}


// Checks the direct solution of a linear system by QRMatrix
template<class MatrixType>
void verification_QRMatrixSolve
(
    MatrixType& A
)
{
    typedef SquareMatrix<typename MatrixType::cmptType> SMatrix;

    // Create working copies of matrix A
    const MatrixType Aorig(A);

    // Direct solution of a linear system by QR decomposition
    #if (0 | RUNALL)
    {
        MatrixType A0(A);
        QRMatrix<MatrixType> QRM
        (
            A0,
            QRMatrix<MatrixType>::outputTypes::FULL_QR
        );

        Info<< nl << "# Solution of A*x = b with A = Q*R:" << nl;
        const scalarField residual(A0.m(), 0);
        const scalarField source(A0.m(), 1);

        const scalarField x(QRM.solve(source));
        const scalarField A0xMinusSource(A0*x - source);
        const scalarField xMinusQRinvSource(x - QRM.inv()*source);
        const SMatrix QRinvA0(QRM.inv()*A0);
        const SMatrix identity(A0.m(), I);

        forAll(residual, i)
        {
            isEqual(A0xMinusSource[i], residual[i], verbose);
            isEqual(xMinusQRinvSource[i], residual[i], verbose);
        }
        equal(QRinvA0, identity, verbose);
    }
    #endif
}


// Checks the parallel tall-skinny QR decomposition
template<class Type>
void verification_tsqr
(
    Random& rndGen
)
{
    typedef RectangularMatrix<Type> RMatrix;

    // Size of the full matrix and its partitions
    const label nColsSum = rndGen.position<label>(1, 100);
    const label qParts = rndGen.position<label>(10, 30);
    List<label> mRowsList(qParts, Zero);

    label mRowsSum = 0;
    for (label i = 0; i < qParts; ++i)
    {
        const label mRows = rndGen.position<label>(nColsSum, 10*nColsSum);
        mRowsList[i] = mRows;
        mRowsSum += mRows;
    }

    # if 0
    Info<< "mRowsSum = " << mRowsSum << nl;
    Info<< "nColsSum = " << nColsSum << nl;
    Info<< "mRowsList = " << mRowsList << nl;
    #endif

    if (mRowsSum < nColsSum)
    {
        FatalErrorInFunction
            << "Tall skinny QR decomposition cannot be executed for matrices "
            << "with number of columns (nCols) = " << nColsSum << " > "
            << "than number of rows (mRows) = " << mRowsSum << "."
            << "Yet nCols > mRows is allowed for a submatrix."
            << exit(FatalError);
    }

    // Perform the QR decomposition for the full matrix
    const RMatrix A(makeRandomMatrix<RMatrix>({mRowsSum, nColsSum}, rndGen));
    QRMatrix<RMatrix> QRM
    (
        A,
        QRMatrix<RMatrix>::outputTypes::REDUCED_R,
        QRMatrix<RMatrix>::colPivoting::FALSE
    );
    const RMatrix& R = QRM.R();

    RMatrix RMag(R);
    for (auto& val : RMag)
    {
        val = mag(val);
    }

    // Partition the full matrix,
    // and perform the QR decomposition on each partition
    RMatrix subRs({qParts*nColsSum, nColsSum}, Zero);
    label mRowsIndex = 0;
    for (label i = 0; i < qParts; ++i)
    {
        RMatrix subA({mRowsList[i], nColsSum}, Zero);
        subA = A.subMatrix(mRowsIndex, 0, mRowsList[i]);
        mRowsIndex += mRowsList[i];

        QRMatrix<RMatrix> subQRM
        (
            subA,
            QRMatrix<RMatrix>::outputTypes::REDUCED_R,
            QRMatrix<RMatrix>::colPivoting::FALSE
        );
        const RMatrix& subR = subQRM.R();

        // Append subR
        subRs.subMatrix(i*nColsSum, 0, nColsSum) = subR;
    }

    // Perform the tall-skinny QR decomposition
    QRMatrix<RMatrix> TSQR
    (
        subRs,
        QRMatrix<RMatrix>::outputTypes::REDUCED_R,
        QRMatrix<RMatrix>::colPivoting::FALSE
    );
    const RMatrix& TSQRR = TSQR.R();

    RMatrix TSQRRMag(TSQRR);
    for (auto& val : TSQRRMag)
    {
        val = mag(val);
    }

    // Compare magnitude of R, since R is not unique up to the sign
    equal(RMag, TSQRRMag, verbose);

    #if 0
    InfoNumPyFormat(R, "R");
    InfoNumPyFormat(TSQRR, "TSQRR");
    #endif
}


// Checks the back substitution function
template<class Type>
void verification_backSubstitution
(
    const SquareMatrix<Type>& A,
    const RectangularMatrix<Type>& b
)
{
    // Ax = b
    typedef RectangularMatrix<Type> RMatrix;

    QRMatrix<RMatrix> QRNull;
    const RMatrix X(QRNull.backSubstitution(A, b));

    #if 1
    {
        InfoNumPyFormat(A, "A");
        InfoNumPyFormat(b, "b");
        InfoNumPyFormat(X, "x");
    }
    #endif

    #if (0 | RUNALL)
    {
        Info<< nl << "# A*X = b:" << nl;
        const RMatrix AX(A*X);
        equal(AX, b, verbose, 10, 1e-3, 1e-3);
    }
    #endif
}


// Checks the Moore-Penrose pseudo-inverse function
template<class MatrixType>
void verification_pinv
(
    const MatrixType& A
)
{
    Info<< "### Moore-Penrose pseudo-inverse: pinv()" << nl << nl;
    const MatrixType APlus(pinv(A));

    #if 0
    InfoNumPyFormat(A, "A");
    InfoNumPyFormat(APlus, "A^+");
    #endif

    if (0 | RUNALL)
    {
        Info<< nl << "# (A^+)^+ = A:" << nl;
        const MatrixType APlusPlus(pinv(APlus));
        equal(APlusPlus, A, verbose);
    }

    // The four Moore-Penrose conditions
    // w.wiki/5RL (Retrieved: 29-06-19)
    if (0 | RUNALL)
    {
        Info<< nl << "# A*(A^+)*A = A:" << nl;
        const MatrixType AAPlusA((A*APlus)*A);
        equal(AAPlusA, A, verbose);
    }

    if (0 | RUNALL)
    {
        Info<< nl << "# (A^+)*A*(A^+) = A:" << nl;
        const MatrixType APlusAAPlus((APlus*A)*APlus);
        equal(APlusAAPlus, APlus, verbose);
    }

    if (0 | RUNALL)
    {
        Info<< nl << "# (A*(A^+)).T() = A*(A^+):" << nl;
        const MatrixType AAPlus(A*APlus);
        const MatrixType AAPlusT(AAPlus.T());
        equal(AAPlusT, AAPlus, verbose);
    }

    if (0 | RUNALL)
    {
        Info<< nl << "# ((A^+)*A = ((A^+)*A).T():" << nl;
        const MatrixType APlusA(APlus*A);
        const MatrixType APlusAT(APlusA.T());
        equal(APlusA, APlusAT, verbose);
    }
}


// * * * * * * * * * * * * * * * Main Program  * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    typedef SquareMatrix<scalar> SMatrix;
    typedef RectangularMatrix<scalar> RMatrix;
    typedef SquareMatrix<complex> SCMatrix;
    typedef RectangularMatrix<complex> RCMatrix;

    Info<< setprecision(15);
    Random rndGen(1234);
    label numberOfTests = 10;

    Info<< "### QR Decomposition:" << nl << nl;

    // SquareMatrix<scalar>
    #if (0 | RUNALL)
    {
        horizontalLine();

        Info<< "## " << numberOfTests << " SquareMatrix<scalar> A:" << nl;

        for (label i = 0; i < numberOfTests; ++i)
        {
            const label mRows = rndGen.position<label>(1, 50);
            Info<< nl << nl << "# Random A with random mRows = " << mRows << nl;

            SMatrix A(makeRandomMatrix<SMatrix>({mRows, mRows}, rndGen));
            #if (0 | RUNALL)
            {
                const SMatrix B(A);
                verification_pinv(B);
            }
            #endif
            verification_QRMatrix(A);
            verification_QRMatrixSolve(A);

        }

        horizontalLine();
    }
    #endif


    // RectangularMatrix<scalar>
    #if (0 | RUNALL)
    {
        horizontalLine();

        Info<< "## " << numberOfTests << " RectangularMatrix<scalar> A:" << nl;

        for (label i = 0; i < numberOfTests; ++i)
        {
            const label mRows = rndGen.position<label>(1, 50);
            const label nCols = rndGen.position<label>(1, 50);
            Info<< nl << nl << "# Random matrix A with"
                << " random mRows = " << mRows
                << " random nCols = " << nCols << nl;

            RMatrix A(makeRandomMatrix<RMatrix>({mRows, nCols}, rndGen));

            #if (0 | RUNALL)
            {
                const RMatrix B(A);
                verification_pinv(B);
            }
            #endif

            verification_QRMatrix(A);
        }

        horizontalLine();
    }
    #endif


    // SquareMatrix<complex(scalar, 0.0)>
    #if (0 | RUNALL)
    {
        horizontalLine();

        Info<< "## " << numberOfTests << " SquareMatrix<complex, 0> A:" << nl;

        for (label i = 0; i < numberOfTests; ++i)
        {
            const label mRows = rndGen.position<label>(1, 50);
            Info<< nl << nl << "# Random A with random mRows = " << mRows << nl;

            SCMatrix A({mRows, mRows}, complex(0, 0));

            for (auto& val : A)
            {
                val.Re() = rndGen.GaussNormal<scalar>();
            }

            verification_QRMatrix(A);
        }

        horizontalLine();
    }
    #endif


    // SquareMatrix<complex(scalar, scalar)>
    #if (0 | RUNALL)
    {
        horizontalLine();

        Info<< "## " << numberOfTests << " SquareMatrix<complex> A:" << nl;

        for (label i = 0; i < numberOfTests; ++i)
        {
            const label mRows = rndGen.position<label>(1, 50);
            Info<< nl << nl << "# Random A with random mRows = " << mRows << nl;

            SCMatrix A(makeRandomMatrix<SCMatrix>({mRows, mRows}, rndGen));

            #if (0 | RUNALL)
            {
                const SCMatrix B(A);
                verification_pinv(B);
            }
            #endif

            verification_QRMatrix(A);
        }

        horizontalLine();
    }
    #endif


    // RectangularMatrix<complex(scalar, scalar)>
    #if (0 | RUNALL)
    {
        horizontalLine();

        Info<< "## " << numberOfTests << " RectangularMatrix<complex> A:" << nl;

        for (label i = 0; i < numberOfTests; ++i)
        {
            const label mRows = rndGen.position<label>(1, 50);
            const label nCols = rndGen.position<label>(1, 50);
            Info<< nl << nl << "# Random matrix A with"
                << " random mRows = " << mRows
                << " random nCols = " << nCols << nl;

            RCMatrix A(makeRandomMatrix<RCMatrix>({mRows, nCols}, rndGen));

            #if (0 | RUNALL)
            {
                const RCMatrix B(A);
                verification_pinv(B);
            }
            #endif

            verification_QRMatrix(A);
        }

        horizontalLine();
    }
    #endif


    #if (0 | RUNALL)
    {
        horizontalLine();

        Info<< "### Parallel direct tall-skinny QR decomposition:" << nl;
        /*
        Parallel direct tall-skinny QR decomposition:
            Benson, A. R., Gleich, D. F., & Demmel, J. (2013).
            Direct QR factorizations for tall-and-skinny matrices in MapReduce
            architectures.
            2013 IEEE International Conference on Big Data.
            DOI:10.1109/bigdata.2013.6691583

            Demmel, J., Grigori, L., Hoemmen, M., & Langou, J. (2012).
            Communication-optimal parallel and sequential QR and LU
            factorizations.
            SIAM Journal on Scientific Computing, 34(1), A206-A239.
            DOI:10.1137/080731992
        */
        // (Benson et al. 2013, Fig. 5) & (Demmel et al. 2012, Fig. 2)

        Info<< "## " << numberOfTests << " <scalar> tests:" << nl;
        for (label i = 0; i < numberOfTests; ++i)
        {
            verification_tsqr<scalar>(rndGen);
        }

        Info<< nl << "## " << numberOfTests << " <complex> tests:" << nl;
        for (label i = 0; i < numberOfTests; ++i)
        {
            verification_tsqr<complex>(rndGen);
        }

        horizontalLine();
    }
    #endif


    #if (0 | RUNALL)
    {
        Info<< "### Back substitution:" << nl << nl;

        numberOfTests = 100;

        // <scalar>
        #if (0 | RUNALL)
        {
            horizontalLine();

            Info<< "## " << numberOfTests << " <scalar>:" << nl;

            for (label i = 0; i < numberOfTests; ++i)
            {
                const label mRows = rndGen.position<label>(20, 50);
                const label pCols = rndGen.position<label>(20, 50);
                Info<< nl << nl << "# Random square matrix A with"
                    << " random mRows = " << mRows
                    << " random nCols = " << mRows << nl;

                SMatrix A(makeRandomMatrix<SMatrix>({mRows, mRows}, rndGen));

                // Zeroise the lower diagonal,
                // so that A becomes an upper-triangular matrix
                for (label i = 0; i < mRows; ++i)
                {
                    for (label j = 0; j < i; ++j)
                    {
                        A(i, j) = 0.0;
                    }
                }

                // det(triangularA) = prod(diag(triangularA))
                scalar det = 1;
                const List<scalar> diag0(A.diag());
                for (const auto& val : diag0)
                {
                    det *= val;
                }

                if (1e-8 < det)
                {
                    Info<< "# Random matrix b with"
                        << " random mRows = " << mRows
                        << " random nCols = " << pCols << nl;

                    const RMatrix b
                    (
                        makeRandomMatrix<RMatrix>({mRows, pCols}, rndGen)
                    );
                    Info<< "det(A) = " << det << nl;
                    verification_backSubstitution(A, b);
                }
                else
                {
                    Info<< "Random matrix A is a nearly-singular matrix."
                        << " Skipping back-substitution" << nl;
                }
            }

            horizontalLine();
        }
        #endif

        // <complex>
        #if (0 | RUNALL)
        {
            horizontalLine();

            Info<< "## " << numberOfTests << " <complex>:" << nl;

            for (label i = 0; i < numberOfTests; ++i)
            {
                const label mRows = rndGen.position<label>(20, 50);
                const label pCols = rndGen.position<label>(20, 50);
                Info<< nl << nl << "# Random square matrix A with"
                    << " random mRows = " << mRows
                    << " random nCols = " << mRows << nl;

                SCMatrix A(makeRandomMatrix<SCMatrix>({mRows, mRows}, rndGen));

                // Zeroise the lower diagonal,
                // so that A becomes an upper-triangular matrix
                for (label i = 0; i < mRows; ++i)
                {
                    for (label j = 0; j < i; ++j)
                    {
                        A(i, j) = pTraits<complex>::zero;
                    }
                }

                // det(triangularA) = prod(diag(triangularA))
                complex det = pTraits<complex>::one;
                const List<complex> diag0(A.diag());
                for (const auto& val : diag0)
                {
                    det *= val;
                }

                if (1e-8 < mag(det))
                {
                    Info<< "# Random matrix b with"
                        << " random mRows = " << mRows
                        << " random nCols = " << pCols << nl;

                    const RCMatrix b
                    (
                        makeRandomMatrix<RCMatrix>({mRows, pCols}, rndGen)
                    );
                    Info<< "det(A) = " << det << nl;
                    verification_backSubstitution(A, b);
                }
                else
                {
                    Info<< "Random matrix A is a nearly-singular matrix."
                        << " Skipping back-substitution" << nl;
                }
            }

            horizontalLine();
        }
        #endif
    }
    #endif


    Info<< nl << "End" << nl;

    return 0;
}


// ************************************************************************* //
