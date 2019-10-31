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
    Foam::DiagonalMatrix

Description
    A 2D diagonal matrix of objects of type \<Type\>, size (N x N)

SourceFiles
    DiagonalMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef DiagonalMatrix_H
#define DiagonalMatrix_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Form, class Type> class Matrix;

/*---------------------------------------------------------------------------*\
                       Class DiagonalMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class DiagonalMatrix
:
    public List<Type>
{
public:

    // Constructors

        //- Construct null.
        DiagonalMatrix<Type>();

        //- Construct empty from size
        explicit DiagonalMatrix<Type>(const label n);

        //- Construct from size
        //- initializing all elements to the given value
        DiagonalMatrix<Type>(const label n, const zero);

        //- Construct from size and a value
        DiagonalMatrix<Type>(const label n, const Type& val);

        //- Construct from the diagonal of a Matrix
        template<class Form>
        DiagonalMatrix<Type>(const Matrix<Form, Type>& mat);


    // Member Functions

        //- Invert the diagonal matrix and return itself
        DiagonalMatrix<Type>& invert();

};


// Global functions

//- Return the diagonal Matrix inverse
template<class Type>
DiagonalMatrix<Type> inv(const DiagonalMatrix<Type>& mat);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DiagonalMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
