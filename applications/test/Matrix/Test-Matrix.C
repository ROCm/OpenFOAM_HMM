/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "scalarMatrices.H"
#include "LUscalarMatrix.H"
#include "LLTMatrix.H"
#include "QRMatrix.H"
#include "vector.H"
#include "tensor.H"
#include "IFstream.H"

#include <algorithm>

using namespace Foam;

// Copy values into matrix
template<class Form, class Type>
void assignMatrix
(
    Matrix<Form, Type>& mat,
    std::initializer_list<typename Matrix<Form, Type>::cmptType> list
)
{
    const label nargs = list.size();

    if (nargs != mat.size())
    {
        FatalErrorInFunction
            << "Mismatch in matrix dimension ("
            << mat.m() << ", "
            << mat.n() << ") and number of args (" << nargs << ')' << nl
            << exit(FatalError);
     }

    std::copy(list.begin(), list.end(), mat.begin());
}

// Create matrix with values
template<class MatrixType>
MatrixType makeMatrix
(
    const labelPair& dims,
    std::initializer_list<typename MatrixType::cmptType> list
)
{
    MatrixType mat(dims);

    assignMatrix(mat, list);

    return mat;
}


// Create matrix with values
template<class MatrixType, Foam::label nRows, Foam::label nCols>
MatrixType makeMatrix
(
    std::initializer_list<typename MatrixType::cmptType> list
)
{
    MatrixType mat(labelPair(nRows, nCols));

    assignMatrix(mat, list);

    return mat;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    // SquareMatrix
    {
        Info<< nl << "Test SquareMatrix" << nl;

        SquareMatrix<scalar> square1(3);

        assignMatrix
        (
            square1,
            {
                -3.0, 10.0, -4.0,
                2.0, 3.0, 10.0,
                2.0, 6.0, 1.0
            }
        );

        Info<< "matrix: " << square1
            << " begin: " << uintptr_t(square1.cdata()) << nl;

        // Test makeMatrix
        {
            auto square1b
            (
                makeMatrix<scalarSquareMatrix>
                (
                    {3, 3},
                    {
                        -3.0, 10.0, -4.0,
                        2.0, 3.0, 10.0,
                        2.0, 6.0, 1.0
                    }
                )
            );

            Info<< "makeMatrix: " << square1b << nl;

            auto square1c
            (
                makeMatrix<scalarSquareMatrix, 3, 3>
                ({
                    -3.0, 10.0, -4.0,
                    2.0, 3.0, 10.0,
                    2.0, 6.0, 1.0
                })
            );

            Info<< "makeMatrix: " << square1c << nl;
        }

        //Info<< square1 - 2.0*(-square1) << nl;
        Info<< "min:" << min(square1) << " max:" << max(square1) << nl;
        Info<< "min/max: " << minMax(square1) << nl;

        // Steal matrix contents

        List<scalar> stole(square1.release());

        Info<< "matrix: " << square1
            << " begin: " << uintptr_t(square1.cdata()) << nl;

        Info<< "List: " << stole << nl;

        Info<< "min/max: " << minMax(square1) << nl;

        square1 = 100;

        Info<< "matrix: " << square1
            << " begin: " << uintptr_t(square1.cdata()) << nl;


        SquareMatrix<scalar> square2(3, I);

        square1 = square2;

        Info<< nl << "After copy assign from Identity:" << nl << square1 << nl;

        square1 += 1.25*square2;

        Info<< nl << "After +=" << nl << square1 << nl;

        square1 -= square2.T();

        Info<< nl << "After -=" << nl << square1 << nl;

        square1 *= 10;

        Info<< nl << "After *=" << nl << square1 << nl;

        square1 /= 8;

        Info<< nl << "After /=" << nl << square1 << nl;

        SquareMatrix<scalar> square4;
        square4 = square2;

        Info<< nl << square4 << endl;

        SquareMatrix<scalar> square5;
        square4 = square5;
        Info<< nl << square5 << endl;
    }

    // RectangularMatrix
    {
        Info<< nl << "Test RectangularMatrix" << nl;

        RectangularMatrix<scalar> rm1({5, 6}, 3.1);
        rm1(0, 1) = 4.5;

        RectangularMatrix<scalar> rm1b(rm1.block(2, 2, 0, 0));

        Info // << "Full matrix " << rm1 << nl
            << "block = " << rm1b << endl;
    }

    {
        scalarSymmetricSquareMatrix symmMatrix(3, Zero);

        symmMatrix(0, 0) = 4;
        symmMatrix(1, 0) = 12;
        symmMatrix(1, 1) = 37;
        symmMatrix(2, 0) = -16;
        symmMatrix(2, 1) = -43;
        symmMatrix(2, 2) = 98;

        Info<< "Symmetric Square Matrix = " << symmMatrix << endl;

        Info<< "Inverse = " << inv(symmMatrix) << endl;
        Info<< "Determinant = " << det(symmMatrix) << endl;

        scalarSymmetricSquareMatrix symmMatrix2(symmMatrix);
        LUDecompose(symmMatrix2);

        Info<< "Inverse = " << invDecomposed(symmMatrix2) << endl;
        Info<< "Determinant = " << detDecomposed(symmMatrix2) << endl;

        scalarDiagonalMatrix rhs(3, Zero);
        rhs[0] = 1;
        rhs[1] = 2;
        rhs[2] = 3;

        LUsolve(symmMatrix, rhs);

        Info<< "Decomposition = " << symmMatrix << endl;
        Info<< "Solution = " << rhs << endl;
    }


    scalarSquareMatrix squareMatrix(3, Zero);

    squareMatrix(0, 0) = 4;
    squareMatrix(0, 1) = 12;
    squareMatrix(0, 2) = -16;
    squareMatrix(1, 0) = 12;
    squareMatrix(1, 1) = 37;
    squareMatrix(1, 2) = -43;
    squareMatrix(2, 0) = -16;
    squareMatrix(2, 1) = -43;
    squareMatrix(2, 2) = 98;

    Info<< nl << "Square Matrix = " << squareMatrix << endl;

    const scalarField source(3, 1);

    {
        {
            scalarSquareMatrix sm(squareMatrix);
            Info<< "det = " << det(sm) << endl;
        }

        scalarSquareMatrix sm(squareMatrix);
        labelList rhs(3, 0);
        label sign;
        LUDecompose(sm, rhs, sign);

        Info<< "Decomposition = " << sm << endl;
        Info<< "Pivots = " << rhs << endl;
        Info<< "Sign = " << sign << endl;
        Info<< "det = " << detDecomposed(sm, sign) << endl;
    }

    {
        LUscalarMatrix LU(squareMatrix);
        scalarField x(LU.solve(source));
        Info<< "LU solve residual " << (squareMatrix*x - source) << endl;

        scalarSquareMatrix inv(3);
        LU.inv(inv);
        Info<< "LU inv " << inv << endl;
        Info<< "LU inv*squareMatrix " << (inv*squareMatrix) << endl;
    }

    {
        LLTMatrix<scalar> LLT(squareMatrix);
        scalarField x(LLT.solve(source));
        Info<< "LLT solve residual " << (squareMatrix*x - source) << endl;
    }

    {
        QRMatrix<scalarSquareMatrix> QR(squareMatrix);
        scalarField x(QR.solve(source));

        Info<< "QR solve residual "
            << (squareMatrix*x - source) << endl;

        Info<< "QR inverse solve residual "
            << (x - QR.inv()*source) << endl;

        Info<< "QR inv *squareMatrix " << (QR.inv()*squareMatrix) << endl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
