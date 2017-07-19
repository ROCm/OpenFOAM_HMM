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
Foam::autoPtr<typename Foam::LduMatrix<Type, DType, LUType>::preconditioner>
Foam::LduMatrix<Type, DType, LUType>::preconditioner::New
(
    const solver& sol,
    const dictionary& preconditionerDict
)
{
    const word preconditionerName = preconditionerDict.lookup("preconditioner");

    if (sol.matrix().symmetric())
    {
        auto cstrIter =
            symMatrixConstructorTablePtr_->cfind(preconditionerName);

        if (!cstrIter.found())
        {
            FatalIOErrorInFunction
            (
                preconditionerDict
            )   << "Unknown symmetric matrix preconditioner "
                << preconditionerName << endl << endl
                << "Valid symmetric matrix preconditioners are :" << endl
                << symMatrixConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }

        return autoPtr<typename LduMatrix<Type, DType, LUType>::preconditioner>
        (
            cstrIter()
            (
                sol,
                preconditionerDict
            )
        );
    }
    else if (sol.matrix().asymmetric())
    {
        auto cstrIter =
            asymMatrixConstructorTablePtr_->cfind(preconditionerName);

        if (!cstrIter.found())
        {
            FatalIOErrorInFunction
            (
                preconditionerDict
            )   << "Unknown asymmetric matrix preconditioner "
                << preconditionerName << endl << endl
                << "Valid asymmetric matrix preconditioners are :" << endl
                << asymMatrixConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }

        return autoPtr<typename LduMatrix<Type, DType, LUType>::preconditioner>
        (
            cstrIter()
            (
                sol,
                preconditionerDict
            )
        );
    }
    else
    {
        FatalIOErrorInFunction
        (
            preconditionerDict
        )   << "cannot preconditione incomplete matrix, "
               "no diagonal or off-diagonal coefficient"
            << exit(FatalIOError);

        return autoPtr<typename LduMatrix<Type, DType, LUType>::preconditioner>
        (
            nullptr
        );
    }
}


// ************************************************************************* //
