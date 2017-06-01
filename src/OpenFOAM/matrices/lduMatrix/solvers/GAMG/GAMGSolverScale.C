/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "GAMGSolver.H"
#include "vector2D.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::GAMGSolver::scale
(
    scalarField& field,
    scalarField& Acf,
    const lduMatrix& A,
    const FieldField<Field, scalar>& interfaceLevelBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaceLevel,
    const scalarField& source,
    const direction cmpt
) const
{
    A.Amul
    (
        Acf,
        field,
        interfaceLevelBouCoeffs,
        interfaceLevel,
        cmpt
    );


    const label nCells = field.size();
    scalar* __restrict__ fieldPtr = field.begin();
    const scalar* const __restrict__ sourcePtr = source.begin();
    const scalar* const __restrict__ AcfPtr = Acf.begin();


    scalar scalingFactorNum = 0.0;
    scalar scalingFactorDenom = 0.0;

    for (label i=0; i<nCells; i++)
    {
        scalingFactorNum += sourcePtr[i]*fieldPtr[i];
        scalingFactorDenom += AcfPtr[i]*fieldPtr[i];
    }

    vector2D scalingVector(scalingFactorNum, scalingFactorDenom);
    A.mesh().reduce(scalingVector, sumOp<vector2D>());

    const scalar sf = scalingVector.x()/stabilise(scalingVector.y(), VSMALL);

    if (debug >= 2)
    {
        Pout<< sf << " ";
    }

    const scalarField& D = A.diag();
    const scalar* const __restrict__ DPtr = D.begin();

    for (label i=0; i<nCells; i++)
    {
        fieldPtr[i] = sf*fieldPtr[i] + (sourcePtr[i] - sf*AcfPtr[i])/DPtr[i];
    }
}


// ************************************************************************* //
