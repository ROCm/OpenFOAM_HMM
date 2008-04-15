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

#include "LduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::autoPtr<typename Foam::LduMatrix<Type, DType, LUType>::smoother>
Foam::LduMatrix<Type, DType, LUType>::smoother::New
(
    const word& fieldName,
    const LduMatrix<Type, DType, LUType>& matrix,
    const dictionary& smootherDict
)
{
    word smootherName = smootherDict.lookup("smoother");

    if (matrix.symmetric())
    {
        typename symMatrixConstructorTable::iterator constructorIter =
            symMatrixConstructorTablePtr_->find(smootherName);

        if (constructorIter == symMatrixConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "LduMatrix<Type, DType, LUType>::smoother::New", smootherDict
            )   << "Unknown symmetric matrix smoother " << smootherName
                << endl << endl
                << "Valid symmetric matrix smoothers are :" << endl
                << symMatrixConstructorTablePtr_->toc()
                << exit(FatalIOError);
        }

        return autoPtr<typename LduMatrix<Type, DType, LUType>::smoother>
        (
            constructorIter()
            (
                fieldName,
                matrix
            )
        );
    }
    else if (matrix.asymmetric())
    {
        typename asymMatrixConstructorTable::iterator constructorIter =
            asymMatrixConstructorTablePtr_->find(smootherName);

        if (constructorIter == asymMatrixConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "LduMatrix<Type, DType, LUType>::smoother::New", smootherDict
            )   << "Unknown asymmetric matrix smoother " << smootherName
                << endl << endl
                << "Valid asymmetric matrix smoothers are :" << endl
                << asymMatrixConstructorTablePtr_->toc()
                << exit(FatalIOError);
        }

        return autoPtr<typename LduMatrix<Type, DType, LUType>::smoother>
        (
            constructorIter()
            (
                fieldName,
                matrix
            )
        );
    }
    else
    {
        FatalIOErrorIn
        (
            "LduMatrix<Type, DType, LUType>::smoother::New", smootherDict
        )   << "cannot solve incomplete matrix, no off-diagonal coefficients"
            << exit(FatalIOError);

        return autoPtr<typename LduMatrix<Type, DType, LUType>::smoother>(NULL);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::LduMatrix<Type, DType, LUType>::smoother::smoother
(
    const word& fieldName,
    const LduMatrix<Type, DType, LUType>& matrix
)
:
    fieldName_(fieldName),
    matrix_(matrix)
{}


// ************************************************************************* //
