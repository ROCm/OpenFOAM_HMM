/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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
#include "FixedList.H"


  #ifndef OMP_UNIFIED_MEMORY_REQUIRED
  #pragma omp requires unified_shared_memory
  #define OMP_UNIFIED_MEMORY_REQUIRED
  #endif


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::GAMGSolver::scale
(
    solveScalarField& field,
    solveScalarField& Acf,
    const lduMatrix& A,
    const FieldField<Field, scalar>& interfaceLevelBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaceLevel,
    const solveScalarField& source,
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
    solveScalar* __restrict__ fieldPtr = field.begin();
    const solveScalar* const __restrict__ sourcePtr = source.begin();
    const solveScalar* const __restrict__ AcfPtr = Acf.begin();


    FixedList<solveScalar, 2> scalingFactor(Zero);

    
    solveScalar scalingFactorNum = 0.0, scalingFactorDenom = 0.0;
   

    #pragma omp target teams distribute parallel for reduction(+:scalingFactorNum, scalingFactorDenom) map(tofrom:scalingFactorNum,scalingFactorDenom) if(target:nCells>20000)
    for (label i=0; i<nCells; i++)
    {
        scalingFactorNum += fieldPtr[i]*sourcePtr[i];
        scalingFactorDenom += fieldPtr[i]*AcfPtr[i];
    }
    scalingFactor[0] = scalingFactorNum;
    scalingFactor[1] = scalingFactorDenom;


    A.mesh().reduce(scalingFactor, sumOp<solveScalar>());

    const solveScalar sf =
    (
        scalingFactor[0]
      / stabilise(scalingFactor[1], pTraits<solveScalar>::vsmall)
    );

    if (debug >= 2)
    {
        Pout<< sf << " ";
    }

    const scalarField& D = A.diag();
    const scalar* const __restrict__ DPtr = D.begin();

      #pragma omp target teams distribute parallel for if(target:nCells>20000)
      for (label i=0; i<nCells; i++)
      {
        fieldPtr[i] = sf*fieldPtr[i] + (sourcePtr[i] - sf*AcfPtr[i])/DPtr[i];
      }
}


// ************************************************************************* //
