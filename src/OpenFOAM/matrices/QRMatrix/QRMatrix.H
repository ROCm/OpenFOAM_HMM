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

Class
    Foam::QRMatrix

Description
    QRMatrix (i.e. QR decomposition, QR factorisation or orthogonal-triangular
    decomposition) decomposes a scalar/complex matrix \c A into the following
    matrix product:

    \verbatim
        A = Q*R,
    \endverbatim

    where
     \c Q is a unitary similarity matrix,
     \c R is an upper triangular matrix.

    Reference:
    \verbatim
        QR decomposition:
            mathworld.wolfram.com/QRDecomposition.html (Retrieved:17-06-19)
            w.wiki/52X (Retrieved: 17-06-19)

        QR decomposition with Householder reflector (tag:M):
            Monahan, J. F. (2011).
            Numerical methods of statistics.
            Cambridge: Cambridge University Press.
            DOI:10.1017/CBO9780511977176

        QR decomposition with column pivoting (tag:QSB):
            Quintana-Ortí, G., Sun, X., & Bischof, C. H. (1998).
            A BLAS-3 version of the QR factorization with column pivoting.
            SIAM Journal on Scientific Computing, 19(5), 1486-1494.
            DOI:10.1137/S1064827595296732

        Moore-Penrose inverse algorithm (tags:KP; KPP):
            Katsikis, V. N., & Pappas, D. (2008).
            Fast computing of the Moore-Penrose inverse matrix.
            Electronic Journal of Linear Algebra, 17(1), 637-650.
            DOI:10.13001/1081-3810.1287

            Katsikis, V. N., Pappas, D., & Petralias, A. (2011).
            An improved method for the computation of
            the Moore–Penrose inverse matrix.
            Applied Mathematics and Computation, 217(23), 9828-9834.
            DOI:10.1016/j.amc.2011.04.080

        Tolerance for the Moore-Penrose inverse algorithm (tag:TA):
            Toutounian, F., & Ataei, A. (2009).
            A new method for computing Moore–Penrose inverse matrices.
            Journal of Computational and applied Mathematics, 228(1), 412-417.
            DOI:10.1016/j.cam.2008.10.008

        Operands of QR decomposition:
            mathworld.wolfram.com/UnitaryMatrix.html (Retrieved:13-06-19)
            mathworld.wolfram.com/UpperTriangularMatrix.html (Retrieved:16-06-19)
            mathworld.wolfram.com/PermutationMatrix.html (Retrieved:20-06-19)

        Back substitution:
            mathworld.wolfram.com/GaussianElimination.html (Retrieved:15-06-19)
    \endverbatim

Usage
    Input types:
     - \c A can be a \c SquareMatrix<Type> or \c RectangularMatrix<Type>

    Output types:
     - \c Q is always of the type of the matrix \c A
     - \c R is always of the type of the matrix \c A

    Options for the output forms of \c QRMatrix (for an (m-by-n) input matrix
    \c A with k = min(m, n)):
     - outputTypes::FULL_R:     computes only \c R                   (m-by-n)
     - outputTypes::FULL_QR:    computes both \c R and \c Q          (m-by-m)
     - outputTypes::REDUCED_R:  computes only reduced \c R           (k-by-n)

    Options where to store \c R:
     - storeMethods::IN_PLACE:        replaces input matrix content with \c R
     - storeMethods::OUT_OF_PLACE:    creates new object of \c R

    Options for the computation of column pivoting:
     - colPivoting::FALSE:            switches off column pivoting
     - colPivoting::TRUE:             switches on column pivoting

    Direct solution of linear systems A x = b is possible by solve() alongside
    the following limitations:
     - \c A         = a scalar square matrix
     - output type  = outputTypes::FULL_QR
     - store method = storeMethods::IN_PLACE

Notes
    - QR decomposition is not unique if \c R is not positive diagonal \c R.
    - The option combination:
      - outputTypes::REDUCED_R
      - storeMethods::IN_PLACE
      will not modify the rows of input matrix \c A after its nth row.
    - Both FULL_R and REDUCED_R QR decompositions execute the same number of
      operations. Yet REDUCED_R QR decomposition returns only the first n rows
      of \c R if m > n for an input m-by-n matrix \c A.
    - For m <= n, FULL_R and REDUCED_R will produce the same matrices.

See also
    Test-QRMatrix.C

SourceFiles
    QRMatrix.C
    QRMatrixI.H

\*---------------------------------------------------------------------------*/

#ifndef QRMatrix_H
#define QRMatrix_H

#include "RectangularMatrix.H"
#include "SquareMatrix.H"
#include "complex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class QRMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class MatrixType>
class QRMatrix
{

public:

    typedef typename MatrixType::cmptType cmptType;
    typedef SquareMatrix<cmptType> SMatrix;
    typedef RectangularMatrix<cmptType> RMatrix;

    //- Options for the output matrix forms of QRMatrix
    enum outputTypes : uint8_t
    {
        FULL_R = 1,         //!< computes only \c R
        FULL_QR = 2,        //!< computes both \c R and \c Q
        REDUCED_R = 3       //!< computes only reduced \c R
    };

    //- Options where to store R
    enum storeMethods : uint8_t
    {
        IN_PLACE = 1,       //!< replaces input matrix content with \c R
        OUT_OF_PLACE = 2    //!< creates new object of \c R
    };

    //- Options for the computation of column pivoting
    enum colPivoting : bool
    {
        FALSE = false,      //!< switches off column pivoting
        TRUE  = true        //!< switches on column pivoting
    };


private:

    // Private Data

        //- Selected option for the output matrix forms of QRMatrix
        outputTypes outputType_;

        //- Selected option where to store R
        const storeMethods storeMethod_;

        //- Selected option for the computation of column pivoting
        const colPivoting colPivot_;

        //- Unitary similarity matrix
        MatrixType Q_;

        //- Upper triangular matrix
        MatrixType R_;

        //- Permutation list (if column-pivoting is on)
        labelList P_;


    // Private Member Functions

        //- Compute the QR decomposition without the column pivoting
        void qr
        (
            MatrixType& A
        );

        //- Compute the QR decomposition with the column pivoting
        //  (QSB:Fig. 1)
        void qrPivot
        (
            MatrixType& A
        );

        //- Compute output matrices in selected forms
        //- using Householder reflectors
        //  (M:Section 4)
        inline void applyHouseholder
        (
            MatrixType& A,
            const RMatrix& reflector,
            const label k = 0
        );

        //- Solve the linear system with the Field argument x initialized to
        //- the appropriate transformed source (e.g. Q.T()*source)
        //- and return the solution in x
        template<template<typename> class ListContainer>
        void solvex(ListContainer<cmptType>& x) const;

        //- Solve the linear system with the given source
        //- and return the solution in x
        template<template<typename> class ListContainer>
        void solveImpl
        (
            List<cmptType>& x,
            const ListContainer<cmptType>& source
        ) const;


protected:

    // Protected Member Functions

        //- Compute Householder reflector on a given matrix column, u
        //  (M:Eq.6.814)
        inline RMatrix householderReflector
        (
            RMatrix u
        );

        //- Apply (in-place) Householder reflectors from the left side: u*A
        inline void applyLeftReflector
        (
            MatrixType& A,
            const RMatrix& reflector,
            const label k = 0,
            const label k1 = 0
        );

        //- Apply (in-place) Householder reflectors from the right side: (u*A)*u
        inline void applyRightReflector
        (
            MatrixType& A,
            const RMatrix& reflector,
            const label k = 0
        );


public:

    // Constructors

        //- Construct null
        QRMatrix();

        //- Construct QRMatrix without performing the decomposition
        QRMatrix
        (
            const outputTypes outputType,
            const storeMethods storeMethod = storeMethods::OUT_OF_PLACE,
            const colPivoting colPivot = colPivoting::FALSE
        );

        //- Construct QRMatrix and perform the QR decomposition
        QRMatrix
        (
            MatrixType& A,
            const outputTypes outputType,
            const storeMethods storeMethod = storeMethods::OUT_OF_PLACE,
            const colPivoting colPivot = colPivoting::FALSE
        );

        //- Construct QRMatrix and perform the QR decomposition
        QRMatrix
        (
            const MatrixType& A,
            const outputTypes outputType,
            const colPivoting colPivot = colPivoting::FALSE
        );


    // Member Functions

        // Information

        //- Return the unitary similarity matrix
        //  Includes implicit round-to-zero as mutable operation
        inline const MatrixType& Q() const;

        //- Return the upper triangular matrix
        inline const MatrixType& R() const;

        //- Return the permutation order (P) as a list
        inline const labelList& orderP() const;

        //- Create and return the permutation matrix
        inline SMatrix P() const;


        // Algorithm

        //- Compute QR decomposition according to constructor settings
        void decompose
        (
            MatrixType& A
        );

        //- Compute QR decomposition according to constructor settings
        void decompose
        (
            const MatrixType& A
        );

        //- Solve the linear system with the given source
        //- and return the solution in the argument x
        void solve
        (
            List<cmptType>& x,
            const UList<cmptType>& source
        ) const;

        //- Solve the linear system with the given source
        //- and return the solution in the argument x
        template<class Addr>
        void solve
        (
            List<cmptType>& x,
            const IndirectListBase<cmptType, Addr>& source
        ) const;

        //- Solve the linear system with the given source
        //- and return the solution
        tmp<Field<cmptType>> solve
        (
            const UList<cmptType>& source
        ) const;

        //- Solve the linear system with the given source
        //- and return the solution
        template<class Addr>
        tmp<Field<cmptType>> solve
        (
            const IndirectListBase<cmptType, Addr>& source
        ) const;

        //- Return the inverse of (Q*R), so that solving x = (Q*R).inv()*source
        SMatrix inv() const;

        //- Solve a row-echelon-form linear system starting from the bottom
        //  Back substitution: Ax = b where:
        //  A = Non-singular upper-triangular square matrix (m-by-m)
        //  x = Solution (m-by-p)
        //  b = Source (m-by-p)
        RMatrix backSubstitution
        (
            const SMatrix& A,
            const RMatrix& b
        );
};


// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

//- (Out-of-place) Moore-Penrose inverse of singular/non-singular
//- square/rectangular scalar/complex matrices
//  (KPP:p. 9834; KP:p. 648)
//  The tolerance (i.e. tol) to ensure the R1 matrix full-rank is given as 1e-13
//  in the original paper (KPP:p. 9834).  Nonetheless, the tolerance = 1e-5
//  given by (TA; mentioned in (KPP:p. 9832)) was observed to be more robust
//  for the tested scenarios.
template<class MatrixType>
MatrixType pinv
(
    const MatrixType& A,
    const scalar tol = 1e-5
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "QRMatrixI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "QRMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
