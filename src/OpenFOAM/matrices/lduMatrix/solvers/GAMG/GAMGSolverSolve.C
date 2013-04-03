/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "GAMGSolver.H"
#include "ICCG.H"
#include "BICCG.H"
#include "SubField.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::solverPerformance Foam::GAMGSolver::solve
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt
) const
{
    // Setup class containing solver performance data
    solverPerformance solverPerf(typeName, fieldName_);

    // Calculate A.psi used to calculate the initial residual
    scalarField Apsi(psi.size());
    matrix_.Amul(Apsi, psi, interfaceBouCoeffs_, interfaces_, cmpt);

    // Create the storage for the finestCorrection which may be used as a
    // temporary in normFactor
    scalarField finestCorrection(psi.size());

    // Calculate normalisation factor
    scalar normFactor = this->normFactor(psi, source, Apsi, finestCorrection);

    if (debug >= 2)
    {
        Pout<< "   Normalisation factor = " << normFactor << endl;
    }

    // Calculate initial finest-grid residual field
    scalarField finestResidual(source - Apsi);

    // Calculate normalised residual for convergence test
    solverPerf.initialResidual() = gSumMag
    (
        finestResidual,
        matrix().mesh().comm()
    )/normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();


    // Check convergence, solve if not converged
    if (!solverPerf.checkConvergence(tolerance_, relTol_))
    {
        // Create coarse grid correction fields
        PtrList<scalarField> coarseCorrFields;

        // Create coarse grid sources
        PtrList<scalarField> coarseSources;

        // Create the smoothers for all levels
        PtrList<lduMatrix::smoother> smoothers;

        // Initialise the above data structures
        initVcycle(coarseCorrFields, coarseSources, smoothers);

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
                coarseCorrFields,
                coarseSources,
                cmpt
            );

            // Calculate finest level residual field
            matrix_.Amul(Apsi, psi, interfaceBouCoeffs_, interfaces_, cmpt);
            finestResidual = source;
            finestResidual -= Apsi;

            solverPerf.finalResidual() = gSumMag
            (
                finestResidual,
                matrix().mesh().comm()
            )/normFactor;

            if (debug >= 2)
            {
                solverPerf.print(Info(matrix().mesh().comm()));
            }
        } while
        (
            ++solverPerf.nIterations() < maxIter_
         && !(solverPerf.checkConvergence(tolerance_, relTol_))
        );
    }

    return solverPerf;
}


void Foam::GAMGSolver::Vcycle
(
    const PtrList<lduMatrix::smoother>& smoothers,
    scalarField& psi,
    const scalarField& source,
    scalarField& Apsi,
    scalarField& finestCorrection,
    scalarField& finestResidual,
    PtrList<scalarField>& coarseCorrFields,
    PtrList<scalarField>& coarseSources,
    const direction cmpt
) const
{
    //debug = 2;

    const label coarsestLevel = matrixLevels_.size() - 1;

    // Restrict finest grid residual for the next level up.
    agglomeration_.restrictField(coarseSources[0], finestResidual, 0, true);

    if (debug >= 2 && nPreSweeps_)
    {
        Pout<< "Pre-smoothing scaling factors: ";
    }


    // Residual restriction (going to coarser levels)
    for (label leveli = 0; leveli < coarsestLevel; leveli++)
    {
        Pout<< "Restriction for level:" << leveli << endl;

        if (coarseSources.set(leveli + 1))
        {
            // If the optional pre-smoothing sweeps are selected
            // smooth the coarse-grid field for the restriced source
            if (nPreSweeps_)
            {
                coarseCorrFields[leveli] = 0.0;

                smoothers[leveli + 1].smooth
                (
                    coarseCorrFields[leveli],
                    coarseSources[leveli],
                    cmpt,
                    min
                    (
                        nPreSweeps_ +  preSweepsLevelMultiplier_*leveli,
                        maxPreSweeps_
                    )
                );

                scalarField::subField ACf
                (
                    Apsi,
                    coarseCorrFields[leveli].size()
                );

                // Scale coarse-grid correction field
                // but not on the coarsest level because it evaluates to 1
                if (scaleCorrection_ && leveli < coarsestLevel - 1)
                {
                    int comm = matrixLevels_[leveli].mesh().comm();

                    scale
                    (
                        coarseCorrFields[leveli],
                        const_cast<scalarField&>
                        (
                            ACf.operator const scalarField&()
                        ),
                        matrixLevels_[leveli],
                        comm,
                        interfaceLevelsBouCoeffs_[leveli],
                        interfaceLevels_[leveli],
                        coarseSources[leveli],
                        cmpt
                    );
                }

                // Correct the residual with the new solution
                matrixLevels_[leveli].Amul
                (
                    const_cast<scalarField&>
                    (
                        ACf.operator const scalarField&()
                    ),
                    coarseCorrFields[leveli],
                    interfaceLevelsBouCoeffs_[leveli],
                    interfaceLevels_[leveli],
                    cmpt
                );

                coarseSources[leveli] -= ACf;
            }

            // Residual is equal to source
            agglomeration_.restrictField
            (
                coarseSources[leveli + 1],
                coarseSources[leveli],
                leveli + 1,
                true
            );
        }
    }

    if (debug >= 2 && nPreSweeps_)
    {
        Pout<< endl;
    }


    // Solve Coarsest level with either an iterative or direct solver
    if (coarseCorrFields.set(coarsestLevel))
    {
        Pout<< "Coarsest solve for level:" << coarsestLevel << endl;
        solveCoarsestLevel
        (
            coarseCorrFields[coarsestLevel],
            coarseSources[coarsestLevel]
        );
    }

    if (debug >= 2)
    {
        Pout<< "Post-smoothing scaling factors: ";
    }

    // Smoothing and prolongation of the coarse correction fields
    // (going to finer levels)

    scalarField dummyField(0);

    for (label leveli = coarsestLevel - 1; leveli >= 0; leveli--)
    {
        Pout<< "Smoothing and prolongation for level:" << leveli << endl;

        if (coarseCorrFields.set(leveli))
        {
            Pout<< "Prolonging from " << leveli + 1 << " up to "
                << leveli << endl;
            //// Create a field for the pre-smoothed correction field
            //// as a sub-field of the finestCorrection which is not
            //// currently being used
            //scalarField::subField preSmoothedCoarseCorrField
            //(
            //    finestCorrection,
            //    coarseCorrFields[leveli].size()
            //);
            scalarField preSmoothedCoarseCorrField;

            // Only store the preSmoothedCoarseCorrField if pre-smoothing is
            // used
            if (nPreSweeps_)
            {
                //preSmoothedCoarseCorrField.assign(coarseCorrFields[leveli]);
                preSmoothedCoarseCorrField = coarseCorrFields[leveli];
            }

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

            Pout<< "Doing stuff at level " << leveli << endl;

            //// Create A.psi for this coarse level as a sub-field of Apsi
            //scalarField::subField ACf
            //(
            //    Apsi,
            //    coarseCorrFields[leveli].size()
            //);
            //scalarField& ACfRef =
            //    const_cast<scalarField&>(ACf.operator const scalarField&());
            scalarField ACfRef(coarseCorrFields[leveli].size());


            if (interpolateCorrection_)
            {
                Pout<< "doing interpolate." << endl;
                interpolate
                (
                    coarseCorrFields[leveli],
                    ACfRef,
                    matrixLevels_[leveli],
                    interfaceLevelsBouCoeffs_[leveli],
                    interfaceLevels_[leveli],
                    coarseSources[leveli],
                    cmpt
                );
                Pout<< "done interpolate." << endl;
            }

            // Scale coarse-grid correction field
            // but not on the coarsest level because it evaluates to 1
            if (scaleCorrection_ && leveli < coarsestLevel - 1)
            {
                //int comm =
                //(
                //    matrixLevels_.set(leveli+1)
                //  ? matrixLevels_[leveli+1].mesh().comm()
                //  : matrixLevels_[leveli].mesh().comm()
                //);
                int comm = matrixLevels_[leveli].mesh().comm();


                Pout<< "doing scale with comm:" << comm << endl;
                scale
                (
                    coarseCorrFields[leveli],
                    ACfRef,
                    matrixLevels_[leveli],
                    comm,
                    interfaceLevelsBouCoeffs_[leveli],
                    interfaceLevels_[leveli],
                    coarseSources[leveli],
                    cmpt
                );
                Pout<< "done scale with comm:" << comm << endl;
            }

            // Only add the preSmoothedCoarseCorrField if pre-smoothing is
            // used
            if (nPreSweeps_)
            {
                coarseCorrFields[leveli] += preSmoothedCoarseCorrField;
            }

            Pout<< "doing smooth." << endl;
            smoothers[leveli + 1].smooth
            (
                coarseCorrFields[leveli],
                coarseSources[leveli],
                cmpt,
                min
                (
                    nPostSweeps_ + postSweepsLevelMultiplier_*leveli,
                    maxPostSweeps_
                )
            );
            Pout<< "done smooth." << endl;
        }
    }

    // Prolong the finest level correction
    Pout<< "Doing Prolong to finest level" << endl;
    agglomeration_.prolongField
    (
        finestCorrection,
        coarseCorrFields[0],
        0,
        true            //false               // no proc agglomeration for now
    );
    Pout<< "Done Prolong to finest level" << endl;

    if (interpolateCorrection_)
    {
        Pout<< "doing interpolate on finest level" << endl;
        interpolate
        (
            finestCorrection,
            Apsi,
            matrix_,
            interfaceBouCoeffs_,
            interfaces_,
            finestResidual,
            cmpt
        );
        Pout<< "done interpolate on finest level" << endl;
    }

    if (scaleCorrection_)
    {
        // Scale the finest level correction
        int comm = matrix_.mesh().comm();
        Pout<< "doing scale on finest level with comm:" << comm << endl;
        scale
        (
            finestCorrection,
            Apsi,
            matrix_,
            comm,
            interfaceBouCoeffs_,
            interfaces_,
            finestResidual,
            cmpt
        );
        Pout<< "done scale on finest level with comm:" << comm << endl;
    }

    forAll(psi, i)
    {
        psi[i] += finestCorrection[i];
    }

    Pout<< "Doing smooth on finest level" << endl;
    smoothers[0].smooth
    (
        psi,
        source,
        cmpt,
        nFinestSweeps_
    );
    Pout<< "Done smooth on finest level" << endl;
}


void Foam::GAMGSolver::initVcycle
(
    PtrList<scalarField>& coarseCorrFields,
    PtrList<scalarField>& coarseSources,
    PtrList<lduMatrix::smoother>& smoothers
) const
{
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

            Pout<< "initVCucle level:" << leveli << " nCoarseCells:"
                << nCoarseCells << endl;

            coarseSources.set(leveli, new scalarField(nCoarseCells));

            //if (!matrixLevels_.set(leveli))
            //{
            //    coarseCorrFields.set(leveli, new scalarField(nCoarseCells));
            //}
        }

        if (matrixLevels_.set(leveli))
        {
            const lduMatrix& mat = matrixLevels_[leveli];

            label nCoarseCells = mat.diag().size();
            Pout<< "initVCucle level:" << leveli << " matrix size:"
                << nCoarseCells << endl;

            coarseCorrFields.set(leveli, new scalarField(nCoarseCells));
            //coarseCorrFields.set(leveli, new scalarField(nCoarseCells));
            //coarseSources.set(leveli, new scalarField(nCoarseCells));

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
}


void Foam::GAMGSolver::solveCoarsestLevel
(
    scalarField& coarsestCorrField,
    const scalarField& coarsestSource
) const
{
    const label coarsestLevel = matrixLevels_.size() - 1;

    Pout<< "solveCoarsestLevel :"
        << " coarsestLevel:" << coarsestLevel << endl;


    label coarseComm = matrixLevels_[coarsestLevel].mesh().comm();
    label oldWarn = UPstream::warnComm;
    UPstream::warnComm = coarseComm;

    if (directSolveCoarsest_)
    {
        coarsestCorrField = coarsestSource;
        coarsestLUMatrixPtr_->solve(coarsestCorrField);
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
    //    label oldWarn = UPstream::warnComm;
    //    UPstream::warnComm = solveComm;
    //
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
    //                coarseSolverPerf = BICCG
    //                (
    //                    "coarsestLevelCorr",
    //                    allMatrix,
    //                    procInterfaceLevelsBouCoeffs_[coarsestLevel],
    //                    procInterfaceLevelsIntCoeffs_[coarsestLevel],
    //                    procInterfaceLevels_[coarsestLevel],
    //                    tolerance_,
    //                    relTol_
    //                ).solve
    //                (
    //                    coarsestCorrField,
    //                    coarsestSource
    //                );
    //            }
    //            else
    //            {
    //                coarseSolverPerf = ICCG
    //                (
    //                    "coarsestLevelCorr",
    //                    allMatrix,
    //                    procInterfaceLevelsBouCoeffs_[coarsestLevel],
    //                    procInterfaceLevelsIntCoeffs_[coarsestLevel],
    //                    procInterfaceLevels_[coarsestLevel],
    //                    tolerance_,
    //                    relTol_
    //                ).solve
    //                (
    //                    coarsestCorrField,
    //                    coarsestSource
    //                );
    //            }
    //        }
    //    }
    //
    //    UPstream::warnComm = oldWarn;
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
    //        coarseSolverPerf.print(Info(coarseComm));
    //    }
    //
    //    Pout<< "procAgglom: coarsestSource   :" << coarsestSource << endl;
    //    Pout<< "procAgglom: coarsestCorrField:" << coarsestCorrField << endl;
    //}
    else
    {
        coarsestCorrField = 0;
        solverPerformance coarseSolverPerf;

        if (matrixLevels_[coarsestLevel].asymmetric())
        {
            coarseSolverPerf = BICCG
            (
                "coarsestLevelCorr",
                matrixLevels_[coarsestLevel],
                interfaceLevelsBouCoeffs_[coarsestLevel],
                interfaceLevelsIntCoeffs_[coarsestLevel],
                interfaceLevels_[coarsestLevel],
                tolerance_,
                relTol_
            ).solve
            (
                coarsestCorrField,
                coarsestSource
            );
        }
        else
        {
            coarseSolverPerf = ICCG
            (
                "coarsestLevelCorr",
                matrixLevels_[coarsestLevel],
                interfaceLevelsBouCoeffs_[coarsestLevel],
                interfaceLevelsIntCoeffs_[coarsestLevel],
                interfaceLevels_[coarsestLevel],
                tolerance_,
                relTol_
            ).solve
            (
                coarsestCorrField,
                coarsestSource
            );
        }

        if (debug >= 2)
        {
            coarseSolverPerf.print(Info(coarseComm));
        }

        Pout<< "GC: coarsestSource   :" << coarsestSource << endl;
        Pout<< "GC: coarsestCorrField:" << coarsestCorrField << endl;
    }

    UPstream::warnComm = oldWarn;
}


// ************************************************************************* //
