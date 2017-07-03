/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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
        auto cstrIter = symMatrixConstructorTablePtr_->cfind(smootherName);

        if (!cstrIter.found())
        {
            FatalIOErrorInFunction(smootherDict)
                << "Unknown symmetric matrix smoother " << smootherName
                << endl << endl
                << "Valid symmetric matrix smoothers are :" << endl
                << symMatrixConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }

        return autoPtr<typename LduMatrix<Type, DType, LUType>::smoother>
        (
            cstrIter()
            (
                fieldName,
                matrix
            )
        );
    }
    else if (matrix.asymmetric())
    {
        auto cstrIter = asymMatrixConstructorTablePtr_->cfind(smootherName);

        if (!cstrIter.found())
        {
            FatalIOErrorInFunction(smootherDict)
                << "Unknown asymmetric matrix smoother " << smootherName
                << endl << endl
                << "Valid asymmetric matrix smoothers are :" << endl
                << asymMatrixConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }

        return autoPtr<typename LduMatrix<Type, DType, LUType>::smoother>
        (
            cstrIter()
            (
                fieldName,
                matrix
            )
        );
    }
    else
    {
        FatalIOErrorInFunction(smootherDict)
            << "cannot solve incomplete matrix, no off-diagonal coefficients"
            << exit(FatalIOError);

        return autoPtr<typename LduMatrix<Type, DType, LUType>::smoother>
        (
            nullptr
        );
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
