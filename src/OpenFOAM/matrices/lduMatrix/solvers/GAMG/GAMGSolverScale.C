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

#ifdef USE_HIP
#include <hip/hip_runtime.h>

__global__ 
static void GAMGSolver_scale_kernel_A( Foam::solveScalar* __restrict__ fieldPtr,  
                      const Foam::solveScalar* const __restrict__ sourcePtr,
                      const Foam::solveScalar* const __restrict__ AcfPtr,
                      Foam::label nCells,
                      Foam::solveScalar* results){

  Foam::label i_start = threadIdx.x+blockIdx.x*blockDim.x;
  Foam::label i_shift = blockDim.x*gridDim.x;

  Foam::solveScalar  scalingFactorNum = 0.0, scalingFactorDenom = 0.0;

  __shared__  Foam::solveScalar s_scalingFactorNum[16];
  __shared__  Foam::solveScalar s_scalingFactorDenom[16]; 
  if (threadIdx.x < 16){
    s_scalingFactorNum[threadIdx.x] = 0.0;
    s_scalingFactorDenom[threadIdx.x] = 0.0;
  }

    for (Foam::label i=i_start; i<nCells; i+=i_shift){
        scalingFactorNum += sourcePtr[i]*fieldPtr[i];
        scalingFactorDenom += AcfPtr[i]*fieldPtr[i];
    }

    //reduce within a warp; assume warp size id 64
    for (int i = 32; i > 0 ; i = i/2){
      scalingFactorNum += __shfl_down(scalingFactorNum,i);
      scalingFactorDenom += __shfl_down(scalingFactorDenom,i);
    }
    //reduce across warps of the same threadblock
    if (threadIdx.x%64 == 0) {
      s_scalingFactorNum[threadIdx.x/64] = scalingFactorNum;
      s_scalingFactorDenom[threadIdx.x/64] = scalingFactorDenom;       
    }
    __syncthreads();
    if (threadIdx.x==0){
        for (int i = 1; i < blockDim.x/64; ++i){ 
          scalingFactorNum   += s_scalingFactorNum[i];
          scalingFactorDenom += s_scalingFactorDenom[i];
        }
        atomicAdd(&results[0],scalingFactorNum);
        atomicAdd(&results[1],scalingFactorDenom);
    }
/*
    //reduce across all warps 
    if (threadIdx.x%64 == 0){
        atomicAdd(&results[0],scalingFactorNum);
        atomicAdd(&results[1],scalingFactorDenom);
    }
*/
}


__global__ 
static void GAMGSolver_scale_kernel_B( Foam::solveScalar* __restrict__ fieldPtr,  
                      const Foam::solveScalar* const __restrict__ sourcePtr,
                      const Foam::solveScalar* const __restrict__ AcfPtr,
                      const Foam::scalar* const __restrict__ DPtr,
                      Foam::scalar sf,
                      Foam::label nCells){

  Foam::label i_start = threadIdx.x+blockIdx.x*blockDim.x;
  Foam::label i_shift = blockDim.x*gridDim.x;
  for (Foam::label i=i_start; i<nCells; i+=i_shift)
        fieldPtr[i] = sf*fieldPtr[i] + (sourcePtr[i] - sf*AcfPtr[i])/DPtr[i];
}

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

    
    #ifdef USE_HIP
      solveScalar * results = new  solveScalar[2];
      results[0] = results[1] = 0.0;

      hipLaunchKernelGGL(HIP_KERNEL_NAME(GAMGSolver_scale_kernel_A), (nCells + 1023)/1024, 1024, 0,0, fieldPtr,  
                      sourcePtr, AcfPtr, nCells, results);
      hipDeviceSynchronize();
      scalingFactorNum = results[0];
      scalingFactorDenom = results[1];

      delete[] results;
    #else
    solveScalar scalingFactorNum = 0.0, scalingFactorDenom = 0.0;
   

    #pragma omp target teams distribute parallel for reduction(+:scalingFactorNum, scalingFactorDenom) map(tofrom:scalingFactorNum,scalingFactorDenom) if(target:nCells>200)
    for (label i=0; i<nCells; i++)
    {
        scalingFactorNum += fieldPtr[i]*sourcePtr[i];
        scalingFactorDenom += fieldPtr[i]*AcfPtr[i];
    }
    scalingFactor[0] = scalingFactorNum;
    scalingFactor[1] = scalingFactorDenom;

    #endif

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

    #ifdef USE_HIP
     hipLaunchKernelGGL(HIP_KERNEL_NAME(GAMGSolver_scale_kernel_B), (nCells + 255)/256, 256, 0,0,
                      fieldPtr, sourcePtr, AcfPtr, DPtr, sf, nCells);
     hipDeviceSynchronize();
    #else
      #pragma omp target teams distribute parallel for if(target:nCells>200)
      for (label i=0; i<nCells; i++)
      {
        fieldPtr[i] = sf*fieldPtr[i] + (sourcePtr[i] - sf*AcfPtr[i])/DPtr[i];
      }
    #endif
}


// ************************************************************************* //
