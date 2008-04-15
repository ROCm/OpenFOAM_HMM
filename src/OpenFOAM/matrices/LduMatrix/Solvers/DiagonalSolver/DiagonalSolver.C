/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "DiagonalSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::DiagonalSolver<Type, DType, LUType>::DiagonalSolver
(
    const word& fieldName,
    const LduMatrix<Type, DType, LUType>& matrix,
    const dictionary& solverDict
)
:
    LduMatrix<Type, DType, LUType>::solver
    (
        fieldName,
        matrix,
        solverDict
    )
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
void Foam::DiagonalSolver<Type, DType, LUType>::read
(
    const dictionary&
)
{}


template<class Type, class DType, class LUType>
typename Foam::LduMatrix<Type, DType, LUType>::solverPerformance
Foam::DiagonalSolver<Type, DType, LUType>::solve
(
    Field<Type>& psi
) const
{
    psi = this->matrix_.source()/this->matrix_.diag();

    return typename LduMatrix<Type, DType, LUType>::solverPerformance
    (
        typeName,
        this->fieldName_,
        pTraits<Type>::zero,
        pTraits<Type>::zero,
        0,
        true,
        false
    );
}


// ************************************************************************* //
