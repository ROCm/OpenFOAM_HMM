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

\*---------------------------------------------------------------------------*/

#include "lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineRunTimeSelectionTable(lduMatrix::smoother, symMatrix);
    defineRunTimeSelectionTable(lduMatrix::smoother, asymMatrix);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::lduMatrix::smoother> Foam::lduMatrix::smoother::New
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    Istream& smootherData
)
{
    word smootherName(smootherData);

    if (matrix.symmetric())
    {
        symMatrixConstructorTable::iterator constructorIter =
            symMatrixConstructorTablePtr_->find(smootherName);

        if (constructorIter == symMatrixConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "lduMatrix::smoother::New", smootherData
            )   << "Unknown symmetric matrix smoother " << smootherName
                << endl << endl
                << "Valid symmetric matrix smoothers are :" << endl
                << symMatrixConstructorTablePtr_->toc()
                << exit(FatalIOError);
        }

        return autoPtr<lduMatrix::smoother>
        (
            constructorIter()
            (
                fieldName,
                matrix,
                interfaceBouCoeffs,
                interfaceIntCoeffs,
                interfaces
            )
        );
    }
    else if (matrix.asymmetric())
    {
        asymMatrixConstructorTable::iterator constructorIter =
            asymMatrixConstructorTablePtr_->find(smootherName);

        if (constructorIter == asymMatrixConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "lduMatrix::smoother::New", smootherData
            )   << "Unknown asymmetric matrix smoother " << smootherName
                << endl << endl
                << "Valid asymmetric matrix smoothers are :" << endl
                << asymMatrixConstructorTablePtr_->toc()
                << exit(FatalIOError);
        }

        return autoPtr<lduMatrix::smoother>
        (
            constructorIter()
            (
                fieldName,
                matrix,
                interfaceBouCoeffs,
                interfaceIntCoeffs,
                interfaces
            )
        );
    }
    else
    {
        FatalIOErrorIn
        (
            "lduMatrix::smoother::New", smootherData
        )   << "cannot solve incomplete matrix, no off-diagonal coefficients"
            << exit(FatalIOError);

        return autoPtr<lduMatrix::smoother>(NULL);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lduMatrix::smoother::smoother
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
)
:
    fieldName_(fieldName),
    matrix_(matrix),
    interfaceBouCoeffs_(interfaceBouCoeffs),
    interfaceIntCoeffs_(interfaceIntCoeffs),
    interfaces_(interfaces)
{}


// ************************************************************************* //
