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

Class
    Foam::Matrix

Description
    A templated (m x n) matrix of objects of \<T\>.
    The layout is (mRows x nCols) - row-major order:

    \verbatim
        (0,0)
          +---> j [nCols]
          |
          |
          v
          i [mRows]
    \endverbatim

SourceFiles
    Matrix.C
    MatrixI.H
    MatrixIO.C

\*---------------------------------------------------------------------------*/

#ifndef Matrix_H
#define Matrix_H

#include "uLabel.H"
#include "Field.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Form, class Type> class Matrix;
template<class MatrixType> class ConstMatrixBlock;
template<class MatrixType> class MatrixBlock;


/*---------------------------------------------------------------------------*\
                           Class Matrix Declaration
\*---------------------------------------------------------------------------*/

template<class Form, class Type>
class Matrix
{
    // Private Data

        //- Number of rows and columns in Matrix.
        label mRows_, nCols_;

        //- Vector of values of type Type
        Type* __restrict__ v_;


    // Private Member Functions

        //- Allocate storage for the contents
        inline void doAlloc();

        //- Multiply matrix with a column vector on the right (A * x)
        template<class ListType>
        tmp<Field<Type>> rightMultiplyImpl(const ListType& colVec) const;

        //- Multiply matrix with a row vector on the left (x * A)
        template<class ListType>
        tmp<Field<Type>> leftMultiplyImpl(const ListType& rowVec) const;


public:

    //- Matrix type
    typedef Matrix<Form, Type> mType;

    //- Component type
    typedef Type cmptType;


    // Static Member Functions

        //- Return a null Matrix
        inline static const Matrix<Form, Type>& null();


    // Iterators

        //- Random access iterator for traversing a Matrix
        typedef Type* iterator;

        //- Random access iterator for traversing a Matrix
        typedef const Type* const_iterator;


    // Constructors

        //- Null constructor.
        inline Matrix();

        //- Construct given number of rows and columns.
        Matrix(const label m, const label n);

        //- Construct with given number of rows and columns
        //- initializing all elements to zero
        Matrix(const label m, const label n, const zero);

        //- Construct with given number of rows and columns
        //- initializing all elements to the given value
        Matrix(const label m, const label n, const Type& val);

        //- Copy construct
        Matrix(const Matrix<Form, Type>& mat);

        //- Move construct
        Matrix(Matrix<Form, Type>&& mat);

        //- Copy constructor from matrix of a different form
        template<class Form2>
        explicit Matrix(const Matrix<Form2, Type>& mat);

        //- Construct from a block of another matrix
        template<class MatrixType>
        Matrix(const ConstMatrixBlock<MatrixType>& Mb);

        //- Construct from a block of another matrix
        template<class MatrixType>
        Matrix(const MatrixBlock<MatrixType>& Mb);

        //- Construct from Istream.
        explicit Matrix(Istream& is);

        //- Clone
        inline autoPtr<mType> clone() const;


    //- Destructor
    ~Matrix();


    // Member Functions

    // Access

        //- Return the number of rows.
        inline label m() const noexcept;

        //- Return the number of columns.
        inline label n() const noexcept;

        //- Return the number of elements in matrix (m*n)
        inline label size() const;

        //- Return true if the matrix is empty (ie, size() is zero)
        inline bool empty() const noexcept;

        //- Return const pointer to the first data element, which can also
        //- be used to address into the matrix contents
        inline const Type* cdata() const noexcept;

        //- Return pointer to the first data element, which can also
        //- be used to address into the matrix contents
        inline Type* data() noexcept;

        //- Return const pointer to data in the specified row.
        //  Subscript checking only with FULLDEBUG
        inline const Type* rowData(const label irow) const;

        //- Return pointer to data in the specified row.
        //  Subscript checking only with FULLDEBUG
        inline Type* rowData(const label irow);

        //- Linear addressing const element access
        //  Subscript checking only with FULLDEBUG
        inline const Type& at(const label idx) const;

        //- Linear addressing element access
        //  Subscript checking only with FULLDEBUG
        inline Type& at(const label idx);



    // Block Access

            inline ConstMatrixBlock<mType> block
            (
                const label m,
                const label n,
                const label mStart,
                const label nStart
            ) const;

            template<class VectorSpace>
            inline ConstMatrixBlock<mType> block
            (
                const label mStart,
                const label nStart
            ) const;

            inline ConstMatrixBlock<mType> col
            (
                const label m,
                const label rowStart
            ) const;

            inline ConstMatrixBlock<mType> col
            (
                const label m,
                const label mStart,
                const label nStart
            ) const;


            inline MatrixBlock<mType> block
            (
                const label m,
                const label n,
                const label mStart,
                const label nStart
            );

            template<class VectorSpace>
            inline MatrixBlock<mType> block
            (
                const label mStart,
                const label nStart
            );

            inline MatrixBlock<mType> col
            (
                const label m,
                const label rowStart
            );

            inline MatrixBlock<mType> col
            (
                const label m,
                const label mStart,
                const label nStart
            );

    // Check

        //- Check index i is within valid range (0 ... m-1).
        inline void checki(const label irow) const;

        //- Check index j is within valid range (0 ... n-1).
        inline void checkj(const label jcol) const;

        //- Check that dimensions are positive, non-zero
        inline void checkSize() const;

        //- True if all entries have identical values, and matrix is non-empty
        inline bool uniform() const;


    // Edit

        //- Clear the Matrix, i.e. set sizes to zero.
        void clear();

        //- Release storage management of the Matrix contents by transferring
        //- management to a List
        List<Type> release();

        //- Swap contents
        void swap(Matrix<Form, Type>& mat);

        //- Transfer the contents of the argument Matrix into this Matrix
        //- and annul the argument Matrix.
        void transfer(Matrix<Form, Type>& mat);

        //- Change the matrix dimensions, preserving the elements
        void resize(const label m, const label n);

        //- Change the matrix dimensions, preserving the elements
        inline void setSize(const label m, const label n);

        //- Resize the matrix without reallocating storage (unsafe)
        inline void shallowResize(const label m, const label n);


    // Operations

        //- Return the transpose of the matrix
        Form T() const;

        //- Multiply matrix with a column vector on the right (A * x)
        inline tmp<Field<Type>> rightMultiply
        (
            const UList<Type>& colVec
        ) const;

        //- Multiply matrix with a column vector on the right (A * x)
        template<class Addr>
        inline tmp<Field<Type>> rightMultiply
        (
            const IndirectListBase<Type, Addr>& colVec
        ) const;

        //- Multiply matrix with a row vector on the left (x * A)
        inline tmp<Field<Type>> leftMultiply
        (
            const UList<Type>& rowVec
        ) const;

        //- Multiply matrix with a row vector on the left (x * A)
        template<class Addr>
        inline tmp<Field<Type>> leftMultiply
        (
            const IndirectListBase<Type, Addr>& rowVec
        ) const;


    // Member Operators

        //- Return const pointer to data in the specified row - rowData().
        //  Subscript checking only with FULLDEBUG
        inline const Type* operator[](const label irow) const;

        //- Return pointer to data in the specified row - rowData().
        //  Subscript checking only with FULLDEBUG
        inline Type* operator[](const label irow);

        //- (i, j) const element access operator
        //  Subscript checking only with FULLDEBUG
        inline const Type& operator()(const label irow, const label jcol) const;

        //- (i, j) element access operator
        //  Subscript checking only with FULLDEBUG
        inline Type& operator()(const label irow, const label jcol);

        //- Copy assignment. Takes linear time.
        void operator=(const Matrix<Form, Type>& mat);

        //- Move assignment
        void operator=(Matrix<Form, Type>&& mat);

        //- Assignment to a block of another matrix
        template<class MatrixType>
        void operator=(const ConstMatrixBlock<MatrixType>& Mb);

        //- Assignment to a block of another matrix
        template<class MatrixType>
        void operator=(const MatrixBlock<MatrixType>& Mb);

        //- Assignment of all elements to zero
        void operator=(const zero);

        //- Assignment of all elements to the given value
        void operator=(const Type& val);

        //- Matrix addition
        void operator+=(const Matrix<Form, Type>& other);

        //- Matrix subtraction
        void operator-=(const Matrix<Form, Type>& other);

        //- Matrix scalar multiplication
        void operator*=(const scalar s);

        //- Matrix scalar division
        void operator/=(const scalar s);


    // Random access iterator (non-const)

        //- Return an iterator to begin traversing a Matrix
        inline iterator begin();

        //- Return an iterator to end traversing a Matrix
        inline iterator end();


    // Random access iterator (const)

        //- Return const_iterator to begin traversing a constant Matrix
        inline const_iterator cbegin() const;

        //- Return const_iterator to end traversing a constant Matrix
        inline const_iterator cend() const;

        //- Return const_iterator to begin traversing a constant Matrix
        inline const_iterator begin() const;

        //- Return const_iterator to end traversing a constant Matrix
        inline const_iterator end() const;


    // IO

        //- Read Matrix from Istream, discarding existing contents.
        bool readMatrix(Istream& is);

        //- Write Matrix, with line-breaks in ASCII when length exceeds
        //- shortLen.
        //  Using '0' suppresses line-breaks entirely.
        Ostream& writeMatrix(Ostream& os, const label shortLen=0) const;


    // Housekeeping

        //- Deprecated(2019-04) raw data pointer, const access
        //  \deprecated(2019-04) - use cdata() method
        const Type* FOAM_DEPRECATED_FOR(2019-04, "cdata() method") v() const
        {
            return v_;
        }

        //- Deprecated(2019-04) raw data pointer, non-const access
        //  \deprecated(2019-04) - use data() method
        Type* FOAM_DEPRECATED_FOR(2019-04, "data() method") v()
        {
            return v_;
        }
};


// IOstream Operators

//- Read Matrix from Istream, discarding contents of existing Matrix.
template<class Form, class Type>
Istream& operator>>(Istream& is, Matrix<Form, Type>& mat);

//- Write Matrix to Ostream, as per Matrix::writeMatrix() with
//- default length, which is given by Detail::ListPolicy::short_length
template<class Form, class Type>
Ostream& operator<<(Ostream& os, const Matrix<Form, Type>& mat);


// Global Functions, Operators

//- Find max value in the matrix
template<class Form, class Type>
const Type& max(const Matrix<Form, Type>& mat);

//- Find min value in the matrix
template<class Form, class Type>
const Type& min(const Matrix<Form, Type>& mat);

//- Find the min/max values of the matrix
template<class Form, class Type>
MinMax<Type> minMax(const Matrix<Form, Type>& mat);

//- Matrix negation
template<class Form, class Type>
Form operator-(const Matrix<Form, Type>& mat);

//- Matrix addition
template<class Form, class Type>
Form operator+
(
    const Matrix<Form, Type>& A,
    const Matrix<Form, Type>& B
);


//- Matrix subtraction
template<class Form, class Type>
Form operator-
(
    const Matrix<Form, Type>& A,
    const Matrix<Form, Type>& B
);


//- Scalar multiplication of a matrix
template<class Form, class Type>
Form operator*
(
    const scalar s,
    const Matrix<Form, Type>& mat
);


//- Scalar multiplication of a matrix
template<class Form, class Type>
Form operator*
(
    const Matrix<Form, Type>& mat,
    const scalar s
);


//- Scalar division of a matrix
template<class Form, class Type>
Form operator/
(
    const Matrix<Form, Type>& mat,
    const scalar s
);


//- Matrix-matrix multiplication
template<class Form1, class Form2, class Type>
typename typeOfInnerProduct<Type, Form1, Form2>::type
operator*
(
    const Matrix<Form1, Type>& A,
    const Matrix<Form2, Type>& B
);


//- Matrix-vector multiplication (A * x), where x is a column vector
template<class Form, class Type>
inline tmp<Field<Type>> operator*
(
    const Matrix<Form, Type>& mat,
    const UList<Type>& x
);

//- Matrix-vector multiplication (A * x), where x is a column vector
template<class Form, class Type, class Addr>
inline tmp<Field<Type>> operator*
(
    const Matrix<Form, Type>& mat,
    const IndirectListBase<Type, Addr>& x
);


//- Vector-matrix multiplication (x * A), where x is a row vector
template<class Form, class Type>
inline tmp<Field<Type>> operator*
(
    const UList<Type>& x,
    const Matrix<Form, Type>& mat
);

//- Vector-matrix multiplication (x * A), where x is a row vector
template<class Form, class Type, class Addr>
inline tmp<Field<Type>> operator*
(
    const IndirectListBase<Type, Addr>& x,
    const Matrix<Form, Type>& mat
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MatrixI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Matrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
