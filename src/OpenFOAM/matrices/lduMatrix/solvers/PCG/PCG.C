/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "PCG.H"
#include "PrecisionAdaptor.H"


#ifdef USE_ROCTX
#include <roctracer/roctx.h>
#endif

//LG using OpenMP offloading and HMM
#ifdef USE_OMP
#include <omp.h>
#endif

#ifndef OMP_UNIFIED_MEMORY_REQUIRED
#pragma omp requires unified_shared_memory
#define OMP_UNIFIED_MEMORY_REQUIRED
#endif


#define USM_PCG 


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PCG, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<PCG>
        addPCGSymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PCG::PCG
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::solver
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::solverPerformance Foam::PCG::scalarSolve
(
    solveScalarField& psi,
    const solveScalarField& source,
    const direction cmpt
) const
{
    #ifdef USE_ROCTX
    roctxRangePush("PCG::scalarSolve");
    #endif

    // --- Setup class containing solver performance data
    solverPerformance solverPerf
    (
        lduMatrix::preconditioner::getName(controlDict_) + typeName,
        fieldName_
    );

    label nCells = psi.size();

    solveScalar* __restrict__ psiPtr = psi.begin();
    
    solveScalarField pA(nCells);
    solveScalar* __restrict__ pAPtr = pA.begin();

    solveScalarField wA(nCells);
    solveScalar* __restrict__ wAPtr = wA.begin();

    solveScalar wArA = solverPerf.great_;
    solveScalar wArAold = wArA;

    #ifdef USE_ROCTX
    roctxRangePush("PCG::Amul");
    #endif


    // --- Calculate A.psi
    matrix_.Amul(wA, psi, interfaceBouCoeffs_, interfaces_, cmpt);
    
    #ifdef USE_ROCTX
    roctxRangePop();
    #endif
    
    
    // --- Calculate initial residual field
    solveScalarField rA(source - wA);
    solveScalar* __restrict__ rAPtr = rA.begin();

    matrix().setResidualField
    (
        ConstPrecisionAdaptor<scalar, solveScalar>(rA)(),
        fieldName_,
        true
    );

    // --- Calculate normalisation factor
    solveScalar normFactor = this->normFactor(psi, source, wA, pA);

    if ((log_ >= 2) || (lduMatrix::debug >= 2))
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }


    #ifdef USE_ROCTX
    roctxRangePush("PCG::gSumMag");
    #endif

    // --- Calculate normalised residual norm
    solverPerf.initialResidual() =
        gSumMag(rA, matrix().mesh().comm())
       /normFactor;

    #ifdef USE_ROCTX
    roctxRangePop();
    #endif


    solverPerf.finalResidual() = solverPerf.initialResidual();

    // --- Check convergence, solve if not converged
    if
    (
        minIter_ > 0
     || !solverPerf.checkConvergence(tolerance_, relTol_, log_)
    )
    {
        // --- Select and construct the preconditioner
        autoPtr<lduMatrix::preconditioner> preconPtr =
            lduMatrix::preconditioner::New
            (
                *this,
                controlDict_
            );

        //printf("num devices=%d\n",omp_get_num_devices());


        // --- Solver iteration
        do
        {
            //printf("LG:  in file %s  line %d, nCells = %d\n",__FILE__, __LINE__, nCells);

            // --- Store previous wArA
            wArAold = wArA;

            #ifdef USE_ROCTX
            roctxRangePush("PCG::precondition");
            #endif
            // --- Precondition residual
            preconPtr->precondition(wA, rA, cmpt);
            #ifdef USE_ROCTX
            roctxRangePop();
            #endif

            #ifdef USE_ROCTX
            roctxRangePush("PCG::gSumProd");
            #endif
            // --- Update search directions:
            wArA = gSumProd(wA, rA, matrix().mesh().comm());
            #ifdef USE_ROCTX
            roctxRangePop();
            #endif

            #ifdef USE_ROCTX
            roctxRangePush("PCG::compute pAPtr");
            #endif

            if (solverPerf.nIterations() == 0)
            {
                  #pragma omp target teams distribute parallel for if(target:nCells>20000) //LG1 AMD
                  for (label cell=0; cell<nCells; cell++)
                  {
                      pAPtr[cell] = wAPtr[cell];
                  }
            }
            else
            {
                  solveScalar beta = wArA/wArAold;
                  #pragma omp target teams distribute parallel for  if(target:nCells>20000) //LG1 AMD
                  for (label cell=0; cell<nCells; cell++)
                  {
                      pAPtr[cell] = wAPtr[cell] + beta*pAPtr[cell];
                  }
            }
            #ifdef USE_ROCTX
            roctxRangePop();
            #endif


            #ifdef USE_ROCTX
            roctxRangePush("PCG::compute Amul");
            #endif

            // --- Update preconditioned residual
            matrix_.Amul(wA, pA, interfaceBouCoeffs_, interfaces_, cmpt);
            
            #ifdef USE_ROCTX
            roctxRangePop();
            #endif

            #ifdef USE_ROCTX
            roctxRangePush("PCG::gSumProd");
            #endif
            
            solveScalar wApA = gSumProd(wA, pA, matrix().mesh().comm());
            
            #ifdef USE_ROCTX
            roctxRangePop();
            #endif

            // --- Test for singularity
            if (solverPerf.checkSingularity(mag(wApA)/normFactor)) break;


            // --- Update solution and residual:

            solveScalar alpha = wArA/wApA;
            
            #ifdef USE_ROCTX
            roctxRangePush("PCG::update psi aA");
            #endif

            #pragma omp target teams distribute parallel for  if(target:nCells>20000) //LG1 AMD
            for (label cell=0; cell<nCells; cell++)
            {
                psiPtr[cell] += alpha*pAPtr[cell];
                rAPtr[cell] -= alpha*wAPtr[cell];
            }

            #ifdef USE_ROCTX
            roctxRangePop();
            #endif


            #ifdef USE_ROCTX
            roctxRangePush("PCG::gSumMag");
            #endif
            solverPerf.finalResidual() =
                gSumMag(rA, matrix().mesh().comm())
               /normFactor;
            #ifdef USE_ROCTX
            roctxRangePop();
            #endif

        } while
        (
            (
              ++solverPerf.nIterations() < maxIter_
            && !solverPerf.checkConvergence(tolerance_, relTol_, log_)
            )
         || solverPerf.nIterations() < minIter_
        );
    }

    matrix().setResidualField
    (
        ConstPrecisionAdaptor<scalar, solveScalar>(rA)(),
        fieldName_,
        false
    );

    //LG1  using roctx marker
    #ifdef USE_ROCTX
    roctxRangePop();
    #endif
    return solverPerf;
}



Foam::solverPerformance Foam::PCG::solve
(
    scalarField& psi_s,
    const scalarField& source,
    const direction cmpt
) const
{
    PrecisionAdaptor<solveScalar, scalar> tpsi(psi_s);
    return scalarSolve
    (
        tpsi.ref(),
        ConstPrecisionAdaptor<solveScalar, scalar>(source)(),
        cmpt
    );
}


// ************************************************************************* //
