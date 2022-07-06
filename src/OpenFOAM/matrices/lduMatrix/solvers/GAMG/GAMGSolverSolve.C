/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
#include "SubField.H"
#include "PrecisionAdaptor.H"

  #ifndef OMP_UNIFIED_MEMORY_REQUIRED
  #pragma omp requires unified_shared_memory
  #define OMP_UNIFIED_MEMORY_REQUIRED
  #endif

#ifdef USE_ROCTX
#include <roctx.h>
#endif


#ifdef USE_HIP 
#include <hip/hip_runtime.h>
__global__
static void  GAMGSolver_Vcycle_kernel_A(Foam::solveScalar* __restrict__ psi, 
                                      const Foam::solveScalar* __restrict__ finestCorrection, Foam::label psi_size){
    Foam::label i_start = threadIdx.x+blockIdx.x*blockDim.x;
    Foam::label i_shift = blockDim.x*gridDim.x;
    for (Foam::label i = i_start; i < psi_size; i+=i_shift)
        psi[i] += finestCorrection[i];
}
#endif


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::solverPerformance Foam::GAMGSolver::solve
(
    scalarField& psi_s,
    const scalarField& source,
    const direction cmpt
) const
{

    printf("in GAMGSolver::solve\n");

    #ifdef USE_ROCTX
    roctxRangePush("GAMGSolver::solve");
    #endif

    PrecisionAdaptor<solveScalar, scalar> tpsi(psi_s);
    solveScalarField& psi = tpsi.ref();

    ConstPrecisionAdaptor<solveScalar, scalar> tsource(source);

    // Setup class containing solver performance data
    solverPerformance solverPerf(typeName, fieldName_);

    // Calculate A.psi used to calculate the initial residual
    solveScalarField Apsi(psi.size());

    #ifdef USE_ROCTX
    roctxRangePush("matrix_.Amul");
    
    #endif
    matrix_.Amul(Apsi, psi, interfaceBouCoeffs_, interfaces_, cmpt);
    #ifdef USE_ROCTX
    roctxRangePop();
    #endif


    // Create the storage for the finestCorrection which may be used as a
    // temporary in normFactor
    solveScalarField finestCorrection(psi.size());

    // Calculate normalisation factor
    solveScalar normFactor =
        this->normFactor(psi, tsource(), Apsi, finestCorrection);

    if ((log_ >= 2) || (debug >= 2))
    {
        Pout<< "   Normalisation factor = " << normFactor << endl;
    }

    // Calculate initial finest-grid residual field
    solveScalarField finestResidual(tsource() - Apsi);

    matrix().setResidualField
    (
        ConstPrecisionAdaptor<scalar, solveScalar>(finestResidual)(),
        fieldName_,
        true
    );

    // Calculate normalised residual for convergence test
    solverPerf.initialResidual() = gSumMag
    (
        finestResidual,
        matrix().mesh().comm()
    )/normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();


    // Check convergence, solve if not converged
    if
    (
        minIter_ > 0
     || !solverPerf.checkConvergence(tolerance_, relTol_, log_)
    )
    {
        // Create coarse grid correction fields
        PtrList<solveScalarField> coarseCorrFields;

        // Create coarse grid sources
        PtrList<solveScalarField> coarseSources;

        // Create the smoothers for all levels
        PtrList<lduMatrix::smoother> smoothers;

        // Scratch fields if processor-agglomerated coarse level meshes
        // are bigger than original. Usually not needed
        solveScalarField scratch1;
        solveScalarField scratch2;

        // Initialise the above data structures
        initVcycle
        (
            coarseCorrFields,
            coarseSources,
            smoothers,
            scratch1,
            scratch2
        );

        do
        {
            Vcycle
            (
                smoothers,
                psi,
                source,
                Apsi,
                finestCorrection,
                finestResidual,

                (scratch1.size() ? scratch1 : Apsi),
                (scratch2.size() ? scratch2 : finestCorrection),

                coarseCorrFields,
                coarseSources,
                cmpt
            );

            // Calculate finest level residual field
            #ifdef USE_ROCTX
            roctxRangePush("matrix_.Amul");
            #endif
            matrix_.Amul(Apsi, psi, interfaceBouCoeffs_, interfaces_, cmpt);
            #ifdef USE_ROCTX
            roctxRangePop();
            #endif
            
            finestResidual = tsource();
            finestResidual -= Apsi;

            solverPerf.finalResidual() = gSumMag
            (
                finestResidual,
                matrix().mesh().comm()
            )/normFactor;

            if ((log_ >= 2) || (debug >= 2))
            {
                solverPerf.print(Info.masterStream(matrix().mesh().comm()));
            }
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
        ConstPrecisionAdaptor<scalar, solveScalar>(finestResidual)(),
        fieldName_,
        false
    );

    #ifdef USE_ROCTX
    roctxRangePop();
    #endif

    return solverPerf;
}


void Foam::GAMGSolver::Vcycle
(
    const PtrList<lduMatrix::smoother>& smoothers,
    solveScalarField& psi,
    const scalarField& source,
    solveScalarField& Apsi,
    solveScalarField& finestCorrection,
    solveScalarField& finestResidual,

    solveScalarField& scratch1,
    solveScalarField& scratch2,

    PtrList<solveScalarField>& coarseCorrFields,
    PtrList<solveScalarField>& coarseSources,
    const direction cmpt
) const
{

    printf("in GAMGSolver::Vcycle\n");

    #ifdef USE_ROCTX
    roctxRangePush("GAMGSolver::Vcycle");
    #endif

    //debug = 2;

    const label coarsestLevel = matrixLevels_.size() - 1;

    // Restrict finest grid residual for the next level up.
    agglomeration_.restrictField(coarseSources[0], finestResidual, 0, true);

    if (nPreSweeps_ && ((log_ >= 2) || (debug >= 2)))
    {
        Pout<< "Pre-smoothing scaling factors: ";
    }


    // Residual restriction (going to coarser levels)
    for (label leveli = 0; leveli < coarsestLevel; leveli++)
    {
        if (coarseSources.set(leveli + 1))
        {
            // If the optional pre-smoothing sweeps are selected
            // smooth the coarse-grid field for the restricted source
            if (nPreSweeps_)
            {
                coarseCorrFields[leveli] = 0.0;

                #ifdef USE_ROCTX
                roctxRangePush("GAMGSolver::Vcycle:smoother_scalarSmooth");
                #endif
                smoothers[leveli + 1].scalarSmooth
                (
                    coarseCorrFields[leveli],
                    coarseSources[leveli],  //coarseSource,
                    cmpt,
                    Foam::min
                    (
                        nPreSweeps_ +  preSweepsLevelMultiplier_*leveli,
                        maxPreSweeps_
                    )
                );
                #ifdef USE_ROCTX
                roctxRangePop();
                #endif
                solveScalarField::subField ACf
                (
                    scratch1,
                    coarseCorrFields[leveli].size()
                );

                // Scale coarse-grid correction field
                // but not on the coarsest level because it evaluates to 1
                if (scaleCorrection_ && leveli < coarsestLevel - 1)
                {
                    #ifdef USE_ROCTX
                    roctxRangePush("GAMGSolver::Vcycle:scale");
                    #endif
                    scale
                    (
                        coarseCorrFields[leveli],
                        const_cast<solveScalarField&>
                        (
                            ACf.operator const solveScalarField&()
                        ),
                        matrixLevels_[leveli],
                        interfaceLevelsBouCoeffs_[leveli],
                        interfaceLevels_[leveli],
                        coarseSources[leveli],
                        cmpt
                    );
                    #ifdef USE_ROCTX
                    roctxRangePop();
                    #endif
                }

                // Correct the residual with the new solution
                #ifdef USE_ROCTX
                roctxRangePush("GAMGSolver::Vcycle:mat_Levels_Amul");
                #endif
                matrixLevels_[leveli].Amul
                (
                    const_cast<solveScalarField&>
                    (
                        ACf.operator const solveScalarField&()
                    ),
                    coarseCorrFields[leveli],
                    interfaceLevelsBouCoeffs_[leveli],
                    interfaceLevels_[leveli],
                    cmpt
                );
                #ifdef USE_ROCTX
                roctxRangePop();
                #endif

                coarseSources[leveli] -= ACf;
            }
            #ifdef USE_ROCTX
            roctxRangePush("GAMGSolver::Vcycle:restrictField");
            #endif
            // Residual is equal to source
            agglomeration_.restrictField
            (
                coarseSources[leveli + 1],
                coarseSources[leveli],
                leveli + 1,
                true
            );
            #ifdef USE_ROCTX
            roctxRangePop();
            #endif
        }
    }

    if (nPreSweeps_ && ((log_ >= 2) || (debug >= 2)))
    {
        Pout<< endl;
    }


    // Solve Coarsest level with either an iterative or direct solver
    if (coarseCorrFields.set(coarsestLevel))
    {
        #ifdef USE_ROCTX
        roctxRangePush("GAMGSolver::Vcycle:solveCoarsestLevel");
        #endif
        solveCoarsestLevel
        (
            coarseCorrFields[coarsestLevel],
            coarseSources[coarsestLevel]
        );
        #ifdef USE_ROCTX
        roctxRangePop();
        #endif   
    }

    if ((log_ >= 2) || (debug >= 2))
    {
        Pout<< "Post-smoothing scaling factors: ";
    }

    // Smoothing and prolongation of the coarse correction fields
    // (going to finer levels)

    solveScalarField dummyField(0);

    for (label leveli = coarsestLevel - 1; leveli >= 0; leveli--)
    {
        if (coarseCorrFields.set(leveli))
        {
            // Create a field for the pre-smoothed correction field
            // as a sub-field of the finestCorrection which is not
            // currently being used
            solveScalarField::subField preSmoothedCoarseCorrField
            (
                scratch2,
                coarseCorrFields[leveli].size()
            );

            // Only store the preSmoothedCoarseCorrField if pre-smoothing is
            // used
            if (nPreSweeps_)
            {
                preSmoothedCoarseCorrField = coarseCorrFields[leveli];
            }


            #ifdef USE_ROCTX
            roctxRangePush("GAMGSolver::Vcycle:prolongField");
            #endif
            agglomeration_.prolongField
            (
                coarseCorrFields[leveli],
                (
                    coarseCorrFields.set(leveli + 1)
                  ? coarseCorrFields[leveli + 1]
                  : dummyField              // dummy value
                ),
                leveli + 1,
                true
            );
            #ifdef USE_ROCTX
            roctxRangePop();
            #endif

            // Create A.psi for this coarse level as a sub-field of Apsi
            solveScalarField::subField ACf
            (
               scratch1,
                coarseCorrFields[leveli].size()
            );
            solveScalarField& ACfRef =
                const_cast
                <
                    solveScalarField&
                >(ACf.operator const solveScalarField&());

            if (interpolateCorrection_) //&& leveli < coarsestLevel - 2)
            {

                #ifdef USE_ROCTX
                roctxRangePush("GAMGSolver::Vcycle:interpolate");
                #endif

                if (coarseCorrFields.set(leveli+1))
                {
                    interpolate
                    (
                        coarseCorrFields[leveli],
                        ACfRef,
                        matrixLevels_[leveli],
                        interfaceLevelsBouCoeffs_[leveli],
                        interfaceLevels_[leveli],
                        agglomeration_.restrictAddressing(leveli + 1),
                        coarseCorrFields[leveli + 1],
                        cmpt
                    );
                }
                else
                {
                    interpolate
                    (
                        coarseCorrFields[leveli],
                        ACfRef,
                        matrixLevels_[leveli],
                        interfaceLevelsBouCoeffs_[leveli],
                        interfaceLevels_[leveli],
                        cmpt
                    );
                }
                #ifdef USE_ROCTX
                roctxRangePop();
                #endif
            }

            // Scale coarse-grid correction field
            // but not on the coarsest level because it evaluates to 1
            #ifdef USE_ROCTX
            roctxRangePush("GAMGSolver::Vcycle:scale-coarse");
            #endif
            if
            (
                scaleCorrection_
             && (interpolateCorrection_ || leveli < coarsestLevel - 1)
            )
            {
                scale
                (
                    coarseCorrFields[leveli],
                    ACfRef,
                    matrixLevels_[leveli],
                    interfaceLevelsBouCoeffs_[leveli],
                    interfaceLevels_[leveli],
                    coarseSources[leveli],
                    cmpt
                );
            }
            #ifdef USE_ROCTX
            roctxRangePop();
            #endif
            // Only add the preSmoothedCoarseCorrField if pre-smoothing is
            // used
            if (nPreSweeps_)
            {
                coarseCorrFields[leveli] += preSmoothedCoarseCorrField;
            }
            #ifdef USE_ROCTX
            roctxRangePush("GAMGSolver::Vcycle:smoother_scalarSmooth");
            #endif
            smoothers[leveli + 1].scalarSmooth
            (
                coarseCorrFields[leveli],
                coarseSources[leveli],  //coarseSource,
                cmpt,
                Foam::min
                (
                    nPostSweeps_ + postSweepsLevelMultiplier_*leveli,
                    maxPostSweeps_
                )
            );
            #ifdef USE_ROCTX
            roctxRangePop();
            #endif
        }
    }

    #ifdef USE_ROCTX
    roctxRangePush("GAMGSolver::Vcycle:prolongField");
    #endif
    // Prolong the finest level correction
    agglomeration_.prolongField
    (
        finestCorrection,
        coarseCorrFields[0],
        0,
        true
    );
    #ifdef USE_ROCTX
    roctxRangePop();
    #endif

    if (interpolateCorrection_)
    {
        #ifdef USE_ROCTX
        roctxRangePush("GAMGSolver::Vcycle:interpolateCorrection");
        #endif        
        interpolate
        (
            finestCorrection,
            Apsi,
            matrix_,
            interfaceBouCoeffs_,
            interfaces_,
            agglomeration_.restrictAddressing(0),
            coarseCorrFields[0],
            cmpt
        );
        #ifdef USE_ROCTX
        roctxRangePop();
        #endif
    }

    if (scaleCorrection_)
    {
        #ifdef USE_ROCTX
        roctxRangePush("GAMGSolver::Vcycle:scaleCorrection");
        #endif
        // Scale the finest level correction
        scale
        (
            finestCorrection,
            Apsi,
            matrix_,
            interfaceBouCoeffs_,
            interfaces_,
            finestResidual,
            cmpt
        );
        #ifdef USE_ROCTX
        roctxRangePop();
        #endif
    }

    #ifdef USE_HIP   //LG2  PROBLEM HERE

       solveScalar* __restrict__ psiPtr = psi.begin();
       const solveScalar* __restrict__ finestCorrectionPtr = finestCorrection.begin();

        hipLaunchKernelGGL(HIP_KERNEL_NAME(GAMGSolver_Vcycle_kernel_A),(psi.size() + 255)/256, 256, 0,0,
                           psiPtr, finestCorrectionPtr, (label) psi.size() );
        hipDeviceSynchronize();
    #else
      #pragma omp target teams distribute parallel for if(target:psi.size()>100)
      for (label i = 0; i < psi.size(); ++i)
      //forAll(psi, i)
      {
        psi[i] += finestCorrection[i];
      }
    #endif


    #ifdef USE_ROCTX
    roctxRangePush("GAMGSolver::Vcycle:smoother_scalarSmooth");
    #endif
    smoothers[0].smooth
    (
        psi,
        source,
        cmpt,
        nFinestSweeps_
    );
    #ifdef USE_ROCTX
    roctxRangePop();
    #endif


    #ifdef USE_ROCTX
    roctxRangePop();
    #endif

}


void Foam::GAMGSolver::initVcycle
(
    PtrList<solveScalarField>& coarseCorrFields,
    PtrList<solveScalarField>& coarseSources,
    PtrList<lduMatrix::smoother>& smoothers,
    solveScalarField& scratch1,
    solveScalarField& scratch2
) const
{
    label maxSize = matrix_.diag().size();

    coarseCorrFields.setSize(matrixLevels_.size());
    coarseSources.setSize(matrixLevels_.size());
    smoothers.setSize(matrixLevels_.size() + 1);

    // Create the smoother for the finest level
    smoothers.set
    (
        0,
        lduMatrix::smoother::New
        (
            fieldName_,
            matrix_,
            interfaceBouCoeffs_,
            interfaceIntCoeffs_,
            interfaces_,
            controlDict_
        )
    );

    forAll(matrixLevels_, leveli)
    {
        if (agglomeration_.nCells(leveli) >= 0)
        {
            label nCoarseCells = agglomeration_.nCells(leveli);

            coarseSources.set(leveli, new solveScalarField(nCoarseCells));
        }

        if (matrixLevels_.set(leveli))
        {
            const lduMatrix& mat = matrixLevels_[leveli];

            label nCoarseCells = mat.diag().size();

            maxSize = Foam::max(maxSize, nCoarseCells);

            coarseCorrFields.set(leveli, new solveScalarField(nCoarseCells));

            smoothers.set
            (
                leveli + 1,
                lduMatrix::smoother::New
                (
                    fieldName_,
                    matrixLevels_[leveli],
                    interfaceLevelsBouCoeffs_[leveli],
                    interfaceLevelsIntCoeffs_[leveli],
                    interfaceLevels_[leveli],
                    controlDict_
                )
            );
        }
    }

    if (maxSize > matrix_.diag().size())
    {
        // Allocate some scratch storage
        scratch1.setSize(maxSize);
        scratch2.setSize(maxSize);
    }
}


Foam::dictionary Foam::GAMGSolver::PCGsolverDict
(
    const scalar tol,
    const scalar relTol
) const
{
    dictionary dict(IStringStream("solver PCG; preconditioner DIC;")());
    dict.add("tolerance", tol);
    dict.add("relTol", relTol);

    return dict;
}


Foam::dictionary Foam::GAMGSolver::PBiCGStabSolverDict
(
    const scalar tol,
    const scalar relTol
) const
{
    dictionary dict(IStringStream("solver PBiCGStab; preconditioner DILU;")());
    dict.add("tolerance", tol);
    dict.add("relTol", relTol);

    return dict;
}


void Foam::GAMGSolver::solveCoarsestLevel
(
    solveScalarField& coarsestCorrField,
    const solveScalarField& coarsestSource
) const
{
    const label coarsestLevel = matrixLevels_.size() - 1;

    const label coarseComm = matrixLevels_[coarsestLevel].mesh().comm();

    if (directSolveCoarsest_)
    {
        PrecisionAdaptor<scalar, solveScalar> tcorrField(coarsestCorrField);

        coarsestLUMatrixPtr_->solve
        (
            tcorrField.ref(),
            ConstPrecisionAdaptor<scalar, solveScalar>(coarsestSource)()
        );
    }
    //else if
    //(
    //    agglomeration_.processorAgglomerate()
    // && procMatrixLevels_.set(coarsestLevel)
    //)
    //{
    //    //const labelList& agglomProcIDs = agglomeration_.agglomProcIDs
    //    //(
    //    //    coarsestLevel
    //    //);
    //    //
    //    //scalarField allSource;
    //    //
    //    //globalIndex cellOffsets;
    //    //if (Pstream::myProcNo(coarseComm) == agglomProcIDs[0])
    //    //{
    //    //    cellOffsets.offsets() =
    //    //        agglomeration_.cellOffsets(coarsestLevel);
    //    //}
    //    //
    //    //cellOffsets.gather
    //    //(
    //    //    coarseComm,
    //    //    agglomProcIDs,
    //    //    coarsestSource,
    //    //    allSource
    //    //);
    //    //
    //    //scalarField allCorrField;
    //    //solverPerformance coarseSolverPerf;
    //
    //    label solveComm = agglomeration_.procCommunicator(coarsestLevel);
    //
    //    coarsestCorrField = 0;
    //    solverPerformance coarseSolverPerf;
    //
    //    if (Pstream::myProcNo(solveComm) != -1)
    //    {
    //        const lduMatrix& allMatrix = procMatrixLevels_[coarsestLevel];
    //
    //        {
    //            Pout<< "** Master:Solving on comm:" << solveComm
    //                << " with procs:" << UPstream::procID(solveComm) << endl;
    //
    //            if (allMatrix.asymmetric())
    //            {
    //                coarseSolverPerf = PBiCGStab
    //                (
    //                    "coarsestLevelCorr",
    //                    allMatrix,
    //                    procInterfaceLevelsBouCoeffs_[coarsestLevel],
    //                    procInterfaceLevelsIntCoeffs_[coarsestLevel],
    //                    procInterfaceLevels_[coarsestLevel],
    //                    PBiCGStabSolverDict(tolerance_, relTol_)
    //                ).solve
    //                (
    //                    coarsestCorrField,
    //                    coarsestSource
    //                );
    //            }
    //            else
    //            {
    //                coarseSolverPerf = PCG
    //                (
    //                    "coarsestLevelCorr",
    //                    allMatrix,
    //                    procInterfaceLevelsBouCoeffs_[coarsestLevel],
    //                    procInterfaceLevelsIntCoeffs_[coarsestLevel],
    //                    procInterfaceLevels_[coarsestLevel],
    //                    PCGsolverDict(tolerance_, relTol_)
    //                ).solve
    //                (
    //                    coarsestCorrField,
    //                    coarsestSource
    //                );
    //            }
    //        }
    //    }
    //
    //    Pout<< "done master solve." << endl;
    //
    //    //// Scatter to all processors
    //    //coarsestCorrField.setSize(coarsestSource.size());
    //    //cellOffsets.scatter
    //    //(
    //    //    coarseComm,
    //    //    agglomProcIDs,
    //    //    allCorrField,
    //    //    coarsestCorrField
    //    //);
    //
    //    if (debug >= 2)
    //    {
    //        coarseSolverPerf.print(Info.masterStream(coarseComm));
    //    }
    //
    //    Pout<< "procAgglom: coarsestSource   :" << coarsestSource << endl;
    //    Pout<< "procAgglom: coarsestCorrField:" << coarsestCorrField << endl;
    //}
    else
    {
        coarsestCorrField = 0;
        const solverPerformance coarseSolverPerf
        (
            coarsestSolverPtr_->scalarSolve
            (
                coarsestCorrField,
                coarsestSource
            )
        );

        if ((log_ >= 2) || debug)
        {
            coarseSolverPerf.print(Info.masterStream(coarseComm));
        }
    }
}


// ************************************************************************* //
