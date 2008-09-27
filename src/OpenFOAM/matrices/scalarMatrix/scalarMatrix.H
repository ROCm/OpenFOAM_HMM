/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::scalarMatrix

Description
    Foam::scalarMatrix

SourceFiles
    scalarMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef scalarMatrix_H
#define scalarMatrix_H

#include "Matrix.H"
#include "DiagonalMatrix.H"
#include "scalarField.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class scalarMatrix Declaration
\*---------------------------------------------------------------------------*/

class scalarMatrix
:
    public Matrix<scalar>
{

public:

    // Constructors

        //- Construct null
        scalarMatrix();

        //- Construct given size
        scalarMatrix(const label);

        //- Construct from Matrix<scalar>
        scalarMatrix(const Matrix<scalar>&);

        //- Construct from Istream
        scalarMatrix(Istream&);


    // Member Functions

        //- Solve the matrix using Gaussian elimination with pivoting,
        //  returning the solution in the source
        template<class Type>
        static void solve(Matrix<scalar>& matrix, Field<Type>& source);

        //- Solve the matrix using Gaussian elimination with pivoting
        //  and return the solution
        template<class Type>
        void solve(Field<Type>& psi, const Field<Type>& source) const;


        //- LU decompose the matrix with pivoting
        static void LUDecompose
        (
            Matrix<scalar>& matrix,
            labelList& pivotIndices
        );

        //- LU back-substitution with given source, returning the solution
        //  in the source
        template<class Type>
        static void LUBacksubstitute
        (
            const Matrix<scalar>& luMmatrix,
            const labelList& pivotIndices,
            Field<Type>& source
        );

        //- Solve the matrix using LU decomposition with pivoting
        //  returning the LU form of the matrix and the solution in the source
        template<class Type>
        static void LUsolve(Matrix<scalar>& matrix, Field<Type>& source);
};


// Global functions

void multiply
(
    Matrix<scalar>& answer,         // value changed in return
    const Matrix<scalar>& A,
    const Matrix<scalar>& B
);

void multiply
(
    Matrix<scalar>& answer,         // value changed in return
    const Matrix<scalar>& A,
    const Matrix<scalar>& B,
    const Matrix<scalar>& C
);

void multiply
(
    Matrix<scalar>& answer,         // value changed in return
    const Matrix<scalar>& A,
    const DiagonalMatrix<scalar>& B,
    const Matrix<scalar>& C
);

//- Return the inverse of matrix A using SVD
Matrix<scalar> SVDinv(const Matrix<scalar>& A, scalar minCondition = 0);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "scalarMatrixTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
