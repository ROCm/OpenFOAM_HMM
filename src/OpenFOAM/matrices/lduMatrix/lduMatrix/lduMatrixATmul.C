/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
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

Description
    Multiply a given vector (second argument) by the matrix or its transpose
    and return the result in the first argument.

\*---------------------------------------------------------------------------*/

#include "lduMatrix.H"



#ifdef USE_ROCTX
#include <roctx.h>
#endif

#ifdef USE_OMP
#include <omp.h>


  #if defined(WM_SP)
  #define _FP_TYPE_scalar float
  #define _FP_TYPE_solve_scalar float
  #elif defined(WM_SPDP)
  #define _FP_TYPE_scalar float
  #define _FP_TYPE_solve_scalar double
  #elif defined(WM_DP)
  #define _FP_TYPE_scalar double
  #define _FP_TYPE_solve_scalar double
  #endif



  #ifndef OMP_UNIFIED_MEMORY_REQUIRED
  #pragma omp requires unified_shared_memory
  #define OMP_UNIFIED_MEMORY_REQUIRED
  #endif
#endif


#ifdef USE_HIP

  #if defined(WM_SP)
  #define _FP_TYPE_scalar float
  #define _FP_TYPE_solve_scalar float
  #elif defined(WM_SPDP)
  #define _FP_TYPE_scalar float
  #define _FP_TYPE_solve_scalar double
  #elif defined(WM_DP)
  #define _FP_TYPE_scalar double
  #define _FP_TYPE_solve_scalar double
  #endif

#include <hip/hip_runtime.h>
__global__
static void lduMatrixATmul_kernel_A(const Foam::scalar *const __restrict X, const Foam::scalar *const __restrict Y, Foam::solveScalar* Z, Foam::label nCells){
  Foam::label i_start = threadIdx.x+blockIdx.x*blockDim.x;
  Foam::label i_shift = blockDim.x*gridDim.x;

  for (Foam::label cell=i_start; cell<nCells; cell+=i_shift){
      Z[cell] = X[cell]*Y[cell];
  }
}

__global__
static void lduMatrixATmul_kernel_B(const Foam::scalar* const __restrict__ lowerPtr, const Foam::scalar* const __restrict__ upperPtr, 
              const Foam::label* const __restrict__ lPtr, const Foam::label* const __restrict__ uPtr, 
              const Foam::solveScalar* const __restrict__ psiPtr, Foam::solveScalar* __restrict__ ApsiPtr, Foam::label nFaces ){
  
  Foam::label i_start = threadIdx.x+blockIdx.x*blockDim.x;
  Foam::label i_shift = blockDim.x*gridDim.x;

  //forcing to use device scop atomics as those are much faster 
  // then the system scop ;
  // ApsiPtr  must point to coarse-grained memory - hence the special API
  for (Foam::label face=i_start; face<nFaces; face+=i_shift){
      /*atomicAdd*/unsafeAtomicAdd(  &ApsiPtr[uPtr[face]], lowerPtr[face]*psiPtr[lPtr[face]] );
      /*atomicAdd*/unsafeAtomicAdd(  &ApsiPtr[lPtr[face]], upperPtr[face]*psiPtr[uPtr[face]] );  
  }
}

__global__
static void lduMatrixATmul_kernel_C( Foam::scalar * __restrict__ ApsiPtr, const Foam::scalar *const __restrict__ ApsiPtr_work_array, Foam::label nCells){
  Foam::label i_start = threadIdx.x+blockIdx.x*blockDim.x;
  Foam::label i_shift = blockDim.x*gridDim.x;

  for (Foam::label cell=i_start; cell<nCells; cell+=i_shift){
      ApsiPtr[cell] = ApsiPtr_work_array[cell];
  }
}

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::lduMatrix::Amul
(
    solveScalarField& Apsi,
    const tmp<solveScalarField>& tpsi,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    solveScalar* __restrict__ ApsiPtr = Apsi.begin();

    const solveScalarField& psi = tpsi();
    const solveScalar* const __restrict__ psiPtr = psi.begin();

    const scalar* const __restrict__ diagPtr = diag().begin();

    const label* const __restrict__ uPtr = lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = lduAddr().lowerAddr().begin();

    const scalar* const __restrict__ upperPtr = upper().begin();
    const scalar* const __restrict__ lowerPtr = lower().begin();

    const label startRequest = Pstream::nRequests();

    // Initialise the update of interfaced interfaces
    initMatrixInterfaces
    (
        true,
        interfaceBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt
    );


    //printf("LG:  in Amul  file = %s line = %d, ApsiPtr = %p\n",__FILE__,__LINE__, ApsiPtr );
    
    const label nCells = diag().size();
   
    #ifndef USE_HIP

      #ifdef USE_OMP
      //_FP_TYPE_solve_scalar*  ApsiPtr_work_array = (_FP_TYPE_solve_scalar*) omp_target_alloc(sizeof(_FP_TYPE_solve_scalar)*nCells, omp_get_default_device() );
      solveScalar*  ApsiPtr_work_array; 
      if (nCells>0)
        ApsiPtr_work_array  = (solveScalar*) omp_target_alloc(sizeof(solveScalar)*nCells, omp_get_default_device() );
      else
        ApsiPtr_work_array  = (solveScalar*) omp_target_alloc(sizeof(solveScalar)*nCells, omp_get_initial_device() );
      #endif
    #endif

    #ifdef USE_HIP

         solveScalar*  ApsiPtr_work_array; 
         hipMalloc( (void**) &ApsiPtr_work_array, sizeof(solveScalar)*nCells );

        //currently atomics are slow in a fine grain memory. converting to coarsegrained memory so we can use fast atomics
        // bug in applying coarsening to the memory pages already in the GPU
        //hipMemAdvise ( (void*) ApsiPtr, sizeof(solveScalar)*nCells, hipMemAdviseSetCoarseGrain, 0);
    #endif
    //printf("LG:  in Amul  file = %s line = %d\n",__FILE__,__LINE__ );

    
    #ifdef USE_HIP
     hipLaunchKernelGGL(HIP_KERNEL_NAME(lduMatrixATmul_kernel_A), (nCells + 255)/256, 256, 0,0, diagPtr, psiPtr, ApsiPtr_work_array, nCells );
     //hipDeviceSynchronize();
    #else

    #pragma omp target teams distribute parallel for //if(target:nCells>2000)
    for (label cell=0; cell<nCells; cell++)
    {

        //ApsiPtr[cell] = diagPtr[cell]*psiPtr[cell];
        #ifdef USE_OMP
           ApsiPtr_work_array[cell] = diagPtr[cell]*psiPtr[cell];
        #else
          ApsiPtr[cell] = diagPtr[cell]*psiPtr[cell];
        #endif
    }
    #endif
    //printf("LG:  in Amul  file = %s line = %d\n",__FILE__,__LINE__ );

    const label nFaces = upper().size();    

    #ifdef USE_HIP
     hipLaunchKernelGGL(HIP_KERNEL_NAME(lduMatrixATmul_kernel_B), (nCells + 255)/256, 256, 0,0, lowerPtr, upperPtr, 
                                                    lPtr,  uPtr, psiPtr,  ApsiPtr_work_array,  nFaces);
     //hipDeviceSynchronize();
    #else

      #pragma omp target teams distribute parallel for //if(target:nCells>2000) // must be nCells, not nFaces to be consistent 
      for (label face=0; face<nFaces; face++)
      {
        #ifdef USE_OMP
          #pragma omp atomic hint(AMD_fast_fp_atomics)
          ApsiPtr_work_array[uPtr[face]] += lowerPtr[face]*psiPtr[lPtr[face]];
          #pragma omp atomic hint(AMD_fast_fp_atomics)
          ApsiPtr_work_array[lPtr[face]] += upperPtr[face]*psiPtr[uPtr[face]];
        #else 
          #pragma omp atomic 
          ApsiPtr[uPtr[face]] += lowerPtr[face]*psiPtr[lPtr[face]];
          #pragma omp atomic 
          ApsiPtr[lPtr[face]] += upperPtr[face]*psiPtr[uPtr[face]];
        #endif

      }
    #endif


    #ifdef USE_OMP
    #pragma omp target teams distribute parallel for //if(target:nCells>2000)
    for (label cell=0; cell<nCells; cell++)
    {
        ApsiPtr[cell] = ApsiPtr_work_array[cell];
    }

    if (nCells>0)
       omp_target_free(ApsiPtr_work_array, omp_get_default_device() );
    else
       omp_target_free(ApsiPtr_work_array, omp_get_initial_device() );

    #endif
    
    #ifdef USE_HIP
     hipLaunchKernelGGL(HIP_KERNEL_NAME(lduMatrixATmul_kernel_C), (nCells + 255)/256, 256, 0,0, ApsiPtr, ApsiPtr_work_array, nCells );
     hipDeviceSynchronize();
     hipFree(ApsiPtr_work_array);
    #endif


    //printf("LG:  in Amul  file = %s line = %d\n",__FILE__,__LINE__ );

    #ifdef USE_ROCTX
    roctxRangePush("lduMatrix::Amul:updateMatrixInterfaces");
    #endif
    // Update interface interfaces
    updateMatrixInterfaces
    (
        true,
        interfaceBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt,
        startRequest
    );
    #ifdef USE_ROCTX
    roctxRangePop();
    #endif

    tpsi.clear();
}


void Foam::lduMatrix::Tmul
(
    solveScalarField& Tpsi,
    const tmp<solveScalarField>& tpsi,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    solveScalar* __restrict__ TpsiPtr = Tpsi.begin();

    const solveScalarField& psi = tpsi();
    const solveScalar* const __restrict__ psiPtr = psi.begin();

    const scalar* const __restrict__ diagPtr = diag().begin();

    const label* const __restrict__ uPtr = lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = lduAddr().lowerAddr().begin();

    const scalar* const __restrict__ lowerPtr = lower().begin();
    const scalar* const __restrict__ upperPtr = upper().begin();

    const label startRequest = Pstream::nRequests();

    // Initialise the update of interfaced interfaces
    initMatrixInterfaces
    (
        true,
        interfaceIntCoeffs,
        interfaces,
        psi,
        Tpsi,
        cmpt
    );

    const label nCells = diag().size();
    for (label cell=0; cell<nCells; cell++)
    {
        TpsiPtr[cell] = diagPtr[cell]*psiPtr[cell];
    }

    const label nFaces = upper().size();
    for (label face=0; face<nFaces; face++)
    {
        TpsiPtr[uPtr[face]] += upperPtr[face]*psiPtr[lPtr[face]];
        TpsiPtr[lPtr[face]] += lowerPtr[face]*psiPtr[uPtr[face]];
    }

    // Update interface interfaces
    updateMatrixInterfaces
    (
        true,
        interfaceIntCoeffs,
        interfaces,
        psi,
        Tpsi,
        cmpt,
        startRequest
    );

    tpsi.clear();
}


void Foam::lduMatrix::sumA
(
    solveScalarField& sumA,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
) const
{
    solveScalar* __restrict__ sumAPtr = sumA.begin();

    const scalar* __restrict__ diagPtr = diag().begin();

    const label* __restrict__ uPtr = lduAddr().upperAddr().begin();
    const label* __restrict__ lPtr = lduAddr().lowerAddr().begin();

    const scalar* __restrict__ lowerPtr = lower().begin();
    const scalar* __restrict__ upperPtr = upper().begin();

    const label nCells = diag().size();
    const label nFaces = upper().size();

    for (label cell=0; cell<nCells; cell++)
    {
        sumAPtr[cell] = diagPtr[cell];
    }

    for (label face=0; face<nFaces; face++)
    {
        sumAPtr[uPtr[face]] += lowerPtr[face];
        sumAPtr[lPtr[face]] += upperPtr[face];
    }

    // Add the interface internal coefficients to diagonal
    // and the interface boundary coefficients to the sum-off-diagonal
    forAll(interfaces, patchi)
    {
        if (interfaces.set(patchi))
        {
            const labelUList& pa = lduAddr().patchAddr(patchi);
            const scalarField& pCoeffs = interfaceBouCoeffs[patchi];

            forAll(pa, face)
            {
                sumAPtr[pa[face]] -= pCoeffs[face];
            }
        }
    }
}


void Foam::lduMatrix::residual
(
    solveScalarField& rA,
    const solveScalarField& psi,
    const scalarField& source,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    solveScalar* __restrict__ rAPtr = rA.begin();

    const solveScalar* const __restrict__ psiPtr = psi.begin();
    const scalar* const __restrict__ diagPtr = diag().begin();
    const scalar* const __restrict__ sourcePtr = source.begin();

    const label* const __restrict__ uPtr = lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = lduAddr().lowerAddr().begin();

    const scalar* const __restrict__ upperPtr = upper().begin();
    const scalar* const __restrict__ lowerPtr = lower().begin();

    // Parallel boundary initialisation.
    // Note: there is a change of sign in the coupled
    // interface update.  The reason for this is that the
    // internal coefficients are all located at the l.h.s. of
    // the matrix whereas the "implicit" coefficients on the
    // coupled boundaries are all created as if the
    // coefficient contribution is of a source-kind (i.e. they
    // have a sign as if they are on the r.h.s. of the matrix.
    // To compensate for this, it is necessary to turn the
    // sign of the contribution.

    const label startRequest = Pstream::nRequests();

    // Initialise the update of interfaced interfaces
    initMatrixInterfaces
    (
        false,
        interfaceBouCoeffs,
        interfaces,
        psi,
        rA,
        cmpt
    );

    const label nCells = diag().size();
    for (label cell=0; cell<nCells; cell++)
    {
        rAPtr[cell] = sourcePtr[cell] - diagPtr[cell]*psiPtr[cell];
    }


    const label nFaces = upper().size();

    for (label face=0; face<nFaces; face++)
    {
        rAPtr[uPtr[face]] -= lowerPtr[face]*psiPtr[lPtr[face]];
        rAPtr[lPtr[face]] -= upperPtr[face]*psiPtr[uPtr[face]];
    }

    // Update interface interfaces
    updateMatrixInterfaces
    (
        false,
        interfaceBouCoeffs,
        interfaces,
        psi,
        rA,
        cmpt,
        startRequest
    );
}


Foam::tmp<Foam::Field<Foam::solveScalar>> Foam::lduMatrix::residual
(
    const solveScalarField& psi,
    const scalarField& source,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    tmp<solveScalarField> trA(new solveScalarField(psi.size()));
    residual(trA.ref(), psi, source, interfaceBouCoeffs, interfaces, cmpt);
    return trA;
}


Foam::tmp<Foam::scalarField> Foam::lduMatrix::H1() const
{
    auto tH1 = tmp<scalarField>::New(lduAddr().size(), Zero);

    if (lowerPtr_ || upperPtr_)
    {
        scalar* __restrict__ H1Ptr = tH1.ref().begin();

        const label* __restrict__ uPtr = lduAddr().upperAddr().begin();
        const label* __restrict__ lPtr = lduAddr().lowerAddr().begin();

        const scalar* __restrict__ lowerPtr = lower().begin();
        const scalar* __restrict__ upperPtr = upper().begin();

        const label nFaces = upper().size();

        for (label face=0; face<nFaces; face++)
        {
            H1Ptr[uPtr[face]] -= lowerPtr[face];
            H1Ptr[lPtr[face]] -= upperPtr[face];
        }
    }

    return tH1;
}


// ************************************************************************* //
