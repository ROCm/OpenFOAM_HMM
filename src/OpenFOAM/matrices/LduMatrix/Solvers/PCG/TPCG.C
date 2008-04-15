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

#include "TPCG.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::TPCG<Type, DType, LUType>::TPCG
(
    const word& fieldName,
    const LduMatrix<Type, DType, LUType>& matrix,
    const dictionary& solverDict
)
:
    LduMatrix<Type, DType, LUType>::solver
    (
        fieldName,
        matrix,
        solverDict
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
typename Foam::LduMatrix<Type, DType, LUType>::solverPerformance
Foam::TPCG<Type, DType, LUType>::solve(Field<Type>& psi) const
{
    word preconditionerName(this->controlDict_.lookup("preconditioner"));

    // --- Setup class containing solver performance data
    typename LduMatrix<Type, DType, LUType>::solverPerformance solverPerf
    (
        preconditionerName + typeName,
        this->fieldName_
    );

    register label nCells = psi.size();

    Type* __restrict__ psiPtr = psi.begin();

    Field<Type> pA(nCells);
    Type* __restrict__ pAPtr = pA.begin();

    Field<Type> wA(nCells);
    Type* __restrict__ wAPtr = wA.begin();

    Type wArA = this->matrix_.great_*pTraits<Type>::one;
    Type wArAold = wArA;

    // --- Calculate A.psi
    this->matrix_.Amul(wA, psi);

    // --- Calculate initial residual field
    Field<Type> rA(this->matrix_.source() - wA);
    Type* __restrict__ rAPtr = rA.begin();

    // --- Calculate normalisation factor
    Type normFactor = this->normFactor(psi, wA, pA);

    if (LduMatrix<Type, DType, LUType>::debug >= 2)
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

    // --- Calculate normalised residual norm
    solverPerf.initialResidual() = cmptDivide(gSumCmptMag(rA), normFactor);
    solverPerf.finalResidual() = solverPerf.initialResidual();

    // --- Check convergence, solve if not converged
    if (!solverPerf.converged(this->tolerance_, this->relTol_))
    {
        // --- Select and construct the preconditioner
        autoPtr<typename LduMatrix<Type, DType, LUType>::preconditioner>
        preconPtr = LduMatrix<Type, DType, LUType>::preconditioner::New
        (
            *this,
            this->controlDict_
        );

        // --- Solver iteration
        do
        {
            // --- Store previous wArA
            wArAold = wArA;

            // --- Precondition residual
            preconPtr->precondition(wA, rA);

            // --- Update search directions:
            wArA = gSumCmptProd(wA, rA);

            if (solverPerf.nIterations() == 0)
            {
                for (register label cell=0; cell<nCells; cell++)
                {
                    pAPtr[cell] = wAPtr[cell];
                }
            }
            else
            {
                Type beta = cmptDivide
                (
                    wArA,
                    stabilise(wArAold, this->matrix_.vsmall_)
                );

                for (register label cell=0; cell<nCells; cell++)
                {
                    pAPtr[cell] = wAPtr[cell] + cmptMultiply(beta, pAPtr[cell]);
                }
            }


            // --- Update preconditioned residual
            this->matrix_.Amul(wA, pA);

            Type wApA = gSumCmptProd(wA, pA);


            // --- Test for singularity
            if (solverPerf.singular(cmptDivide(cmptMag(wApA), normFactor)))
            {
                break;
            }


            // --- Update solution and residual:

            Type alpha = cmptDivide
            (
                wArA,
                stabilise(wApA, this->matrix_.vsmall_)
            );

            for (register label cell=0; cell<nCells; cell++)
            {
                psiPtr[cell] += cmptMultiply(alpha, pAPtr[cell]);
                rAPtr[cell] -= cmptMultiply(alpha, wAPtr[cell]);
            }

            solverPerf.finalResidual() =
                cmptDivide(gSumCmptMag(rA), normFactor);

        } while
        (
            solverPerf.nIterations()++ < this->maxIter_
        && !(solverPerf.converged(this->tolerance_, this->relTol_))
        );
    }

    return solverPerf;
}


// ************************************************************************* //
