/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::SquareMatrix

Description
    A templated 2D square matrix of objects of \<T\>, where the n x n matrix
    dimension is known and used for subscript bounds checking, etc.

SourceFiles
    SquareMatrixI.H
    SquareMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef SquareMatrix_H
#define SquareMatrix_H

#include "Matrix.H"
#include "Identity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Type> class RectangularMatrix;


/*---------------------------------------------------------------------------*\
                        Class SquareMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class SquareMatrix
:
    public Matrix<SquareMatrix<Type>, Type>
{

public:

    // Constructors

        //- Construct null
        inline SquareMatrix();

        //- Construct for given size (rows == cols)
        inline explicit SquareMatrix(const label n);

        //- Construct for given size (rows == cols)
        //- initializing all elements to zero
        inline SquareMatrix(const label n, const zero);

        //- Construct for given size (rows == cols)
        //- initializing all elements to the given value
        inline SquareMatrix(const label n, const Type& val);

        //- Construct for given size (rows == cols)
        //- initializing to the identity matrix
        template<class AnyType>
        inline SquareMatrix(const label n, const Identity<AnyType>);

        template<class AnyType>
        inline SquareMatrix(const labelPair& dims, const Identity<AnyType>);

        //- Construct given number of rows/columns (checked to be equal)
        //  For constructor consistency in rectangular matrices
        inline explicit SquareMatrix(const labelPair& dims);

        //- Construct given number of rows/columns (checked to be equal)
        //- initializing all elements to zero
        //  For constructor consistency with rectangular matrices
        inline SquareMatrix(const labelPair& dims, const zero);

        //- Construct given number of rows/columns (checked to be equal)
        //- initializing all elements to the given value
        //  For constructor consistency with rectangular matrices
        inline SquareMatrix(const labelPair& dims, const Type& val);

        //- Construct given number of rows/columns (checked to be equal)
        //- initializing all elements to zero
        inline SquareMatrix(const label m, const label n, const zero);

        //- Construct from sub-matrix block
        template<class MatrixType>
        inline SquareMatrix(const ConstMatrixBlock<MatrixType>& mat);

        //- Construct from sub-matrix block
        template<class MatrixType>
        inline SquareMatrix(const MatrixBlock<MatrixType>& mat);

        //- Construct as copy of a RectangularMatrix
        //- which is checked to be square
        inline explicit SquareMatrix(const RectangularMatrix<Type>& mat);

        //- Construct from Istream
        inline explicit SquareMatrix(Istream& is);

        //- Clone
        inline autoPtr<SquareMatrix<Type>> clone() const;


    // Member Functions

    // Edit

        //- Resize the matrix preserving the elements
        inline void resize(const label m);

        //- Resize the matrix preserving the elements (compatibility)
        inline void resize(const label m, const label n);

        //- Resize the matrix preserving the elements
        inline void setSize(const label m);

        //- Resize the matrix without reallocating storage (unsafe)
        inline void shallowResize(const label m);

    // Check

        //- Return true if the square matrix is effectively symmetric/Hermitian
        inline bool symmetric() const;

        //- Return true if the square matrix is reduced tridiagonal
        inline bool tridiagonal() const;


    // Member Operators

        //- Assign all elements to zero
        inline void operator=(const zero);

        //- Assign all elements to value
        inline void operator=(const Type& val);

        //- Set to identity matrix
        template<class AnyType>
        void operator=(const Identity<AnyType>);
};


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Return the LU decomposed SquareMatrix det
template<class Type>
scalar detDecomposed(const SquareMatrix<Type>&, const label sign);

//- Return the SquareMatrix det
template<class Type>
scalar det(const SquareMatrix<Type>&);

//- Return the SquareMatrix det and the LU decomposition in the original matrix
template<class Type>
scalar det(SquareMatrix<Type>&);

template<class Type>
class typeOfInnerProduct<Type, SquareMatrix<Type>, SquareMatrix<Type>>
{
public:

    typedef SquareMatrix<Type> type;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SquareMatrixI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SquareMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
