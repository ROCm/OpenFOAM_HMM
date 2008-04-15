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

#include "lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineRunTimeSelectionTable(lduMatrix::preconditioner, symMatrix);
    defineRunTimeSelectionTable(lduMatrix::preconditioner, asymMatrix);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::lduMatrix::preconditioner>
Foam::lduMatrix::preconditioner::New
(
    const solver& sol,
    Istream& preconditionerData
)
{
    word preconditionerName(preconditionerData);

    if (sol.matrix().symmetric())
    {
        symMatrixConstructorTable::iterator constructorIter =
            symMatrixConstructorTablePtr_->find(preconditionerName);

        if (constructorIter == symMatrixConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "lduMatrix::preconditioner::New(const solver&, Istream&)",
                preconditionerData
            )   << "Unknown symmetric matrix preconditioner "
                << preconditionerName << endl << endl
                << "Valid symmetric matrix preconditioners are :" << endl
                << symMatrixConstructorTablePtr_->toc()
                << exit(FatalIOError);
        }

        return autoPtr<lduMatrix::preconditioner>
        (
            constructorIter()
            (
                sol,
                preconditionerData
            )
        );
    }
    else if (sol.matrix().asymmetric())
    {
        asymMatrixConstructorTable::iterator constructorIter =
            asymMatrixConstructorTablePtr_->find(preconditionerName);

        if (constructorIter == asymMatrixConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "lduMatrix::preconditioner::New(const solver&, Istream&)",
                preconditionerData
            )   << "Unknown asymmetric matrix preconditioner "
                << preconditionerName << endl << endl
                << "Valid asymmetric matrix preconditioners are :" << endl
                << asymMatrixConstructorTablePtr_->toc()
                << exit(FatalIOError);
        }

        return autoPtr<lduMatrix::preconditioner>
        (
            constructorIter()
            (
                sol,
                preconditionerData
            )
        );
    }
    else
    {
        FatalIOErrorIn
        (
            "lduMatrix::preconditioner::New(const solver&, Istream&)",
            preconditionerData
        )   << "cannot preconditione incomplete matrix, "
               "no diagonal or off-diagonal coefficient"
            << exit(FatalIOError);

        return autoPtr<lduMatrix::preconditioner>(NULL);
    }
}


// ************************************************************************* //
