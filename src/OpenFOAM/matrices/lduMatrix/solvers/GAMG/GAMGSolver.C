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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GAMGSolver, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<GAMGSolver>
        addGAMGSolverMatrixConstructorToTable_;

    lduMatrix::solver::addasymMatrixConstructorToTable<GAMGSolver>
        addGAMGAsymSolverMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GAMGSolver::GAMGSolver
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
    ),

    // Default values for all controls
    // which may be overridden by those in controlDict
    cacheAgglomeration_(false),
    nPreSweeps_(0),
    preSweepsLevelMultiplier_(1),
    maxPreSweeps_(4),
    nPostSweeps_(2),
    postSweepsLevelMultiplier_(1),
    maxPostSweeps_(4),
    nFinestSweeps_(2),
    interpolateCorrection_(false),
    scaleCorrection_(matrix.symmetric()),
    directSolveCoarsest_(false),
    processorAgglomerate_(false),
    agglomeration_(GAMGAgglomeration::New(matrix_, controlDict_)),

    matrixLevels_(agglomeration_.size()),
    interfaceLevels_(agglomeration_.size()),
    interfaceLevelsBouCoeffs_(agglomeration_.size()),
    interfaceLevelsIntCoeffs_(agglomeration_.size())
{
    readControls();

    forAll(agglomeration_, fineLevelIndex)
    {
        agglomerateMatrix(fineLevelIndex);
    }

    if (matrixLevels_.size())
    {
        if (directSolveCoarsest_)
        {
            const label coarsestLevel = matrixLevels_.size() - 1;

            Pout<< "GAMGSolver :"
                << " coarsestLevel:" << coarsestLevel << endl;

            const lduMesh& coarsestMesh = matrixLevels_[coarsestLevel].mesh();

            label coarseComm = coarsestMesh.comm();
            label oldWarn = UPstream::warnComm;
            UPstream::warnComm = coarseComm;

            Pout<< "Solve direct on coasestmesh (level=" << coarsestLevel
                << ") using communicator " << coarseComm << endl;

            coarsestLUMatrixPtr_.set
            (
                new LUscalarMatrix
                (
                    matrixLevels_[coarsestLevel],
                    interfaceLevelsBouCoeffs_[coarsestLevel],
                    interfaceLevels_[coarsestLevel]
                )
            );

            UPstream::warnComm = oldWarn;
        }
        else if (processorAgglomerate_)
        {
            // Pick a level to processor agglomerate
            label agglomLevel = matrixLevels_.size() - 1;//1;


            // Get mesh and matrix at this level
            const lduMatrix& levelMatrix = matrixLevels_[agglomLevel];
            const lduMesh& levelMesh = levelMatrix.mesh();


            label levelComm = levelMesh.comm();
            label oldWarn = UPstream::warnComm;
            UPstream::warnComm = levelComm;

            Pout<< "Solve generic on mesh (level=" << agglomLevel
                << ") using communicator " << levelComm << endl;

            // Processor restriction map: per processor the coarse processor
            labelList procAgglomMap(UPstream::nProcs(levelComm));
            // Master processor
            labelList masterProcs;
            // Local processors that agglomerate. agglomProcIDs[0] is in
            // masterProc.
            List<int> agglomProcIDs;

            {
                procAgglomMap[0] = 0;
                procAgglomMap[1] = 0;
                procAgglomMap[2] = 1;
                procAgglomMap[3] = 1;

                // Determine the master processors
                Map<label> agglomToMaster(procAgglomMap.size());

                forAll(procAgglomMap, procI)
                {
                    label coarseI = procAgglomMap[procI];

                    Map<label>::iterator fnd = agglomToMaster.find(coarseI);
                    if (fnd == agglomToMaster.end())
                    {
                        agglomToMaster.insert(coarseI, procI);
                    }
                    else
                    {
                        fnd() = max(fnd(), procI);
                    }
                }

                masterProcs.setSize(agglomToMaster.size());
                forAllConstIter(Map<label>, agglomToMaster, iter)
                {
                    masterProcs[iter.key()] = iter();
                }


                // Collect all the processors in my agglomeration
                label myProcID = Pstream::myProcNo(levelComm);
                label myAgglom = procAgglomMap[myProcID];

                // Get all processors agglomerating to the same coarse processor
                agglomProcIDs = findIndices(procAgglomMap, myAgglom);
                // Make sure the master is the first element.
                label index = findIndex
                (
                    agglomProcIDs,
                    agglomToMaster[myAgglom]
                );
                Swap(agglomProcIDs[0], agglomProcIDs[index]);
            }


            Pout<< "procAgglomMap:" << procAgglomMap << endl;
            Pout<< "agglomProcIDs:" << agglomProcIDs << endl;

            // Allocate a communicator for the processor-agglomerated matrix
            label procAgglomComm = UPstream::allocateCommunicator
            (
                levelComm,
                masterProcs
            );
            Pout<< "** Allocated communicator " << procAgglomComm
                << " for indices " << masterProcs
                << " in processor list " << UPstream::procID(levelComm)
                << endl;



            // Gather matrix and mesh onto agglomProcIDs[0]
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            procAgglomerateMatrix
            (
                // Agglomeration information
                procAgglomMap,
                agglomProcIDs,
                procAgglomComm,

                agglomLevel,  // level (coarse, not fine level!)

                // Resulting matrix
                allMatrixPtr_,
                allInterfaceBouCoeffs_,
                allInterfaceIntCoeffs_,
                allPrimitiveInterfaces_,
                allInterfaces_
            );


            UPstream::warnComm = oldWarn;
        }
    }
    else
    {
        FatalErrorIn
        (
            "GAMGSolver::GAMGSolver"
            "("
            "const word& fieldName,"
            "const lduMatrix& matrix,"
            "const FieldField<Field, scalar>& interfaceBouCoeffs,"
            "const FieldField<Field, scalar>& interfaceIntCoeffs,"
            "const lduInterfaceFieldPtrsList& interfaces,"
            "const dictionary& solverControls"
            ")"
        )   << "No coarse levels created, either matrix too small for GAMG"
               " or nCellsInCoarsestLevel too large.\n"
               "    Either choose another solver of reduce "
               "nCellsInCoarsestLevel."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GAMGSolver::~GAMGSolver()
{
    // Clear the the lists of pointers to the interfaces
    forAll(interfaceLevels_, leveli)
    {
        lduInterfaceFieldPtrsList& curLevel = interfaceLevels_[leveli];

        forAll(curLevel, i)
        {
            if (curLevel.set(i))
            {
                delete curLevel(i);
            }
        }
    }

    if (!cacheAgglomeration_)
    {
        delete &agglomeration_;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::GAMGSolver::readControls()
{
    lduMatrix::solver::readControls();

    // we could also consider supplying defaults here too
    controlDict_.readIfPresent("cacheAgglomeration", cacheAgglomeration_);
    controlDict_.readIfPresent("nPreSweeps", nPreSweeps_);
    controlDict_.readIfPresent
    (
        "preSweepsLevelMultiplier",
        preSweepsLevelMultiplier_
    );
    controlDict_.readIfPresent("maxPreSweeps", maxPreSweeps_);
    controlDict_.readIfPresent("nPostSweeps", nPostSweeps_);
    controlDict_.readIfPresent
    (
        "postSweepsLevelMultiplier",
        postSweepsLevelMultiplier_
    );
    controlDict_.readIfPresent("maxPostSweeps", maxPostSweeps_);
    controlDict_.readIfPresent("nFinestSweeps", nFinestSweeps_);
    controlDict_.readIfPresent("interpolateCorrection", interpolateCorrection_);
    controlDict_.readIfPresent("scaleCorrection", scaleCorrection_);
    controlDict_.readIfPresent("directSolveCoarsest", directSolveCoarsest_);
    controlDict_.readIfPresent("processorAgglomerate", processorAgglomerate_);

    Pout<< "GAMGSolver settings :"
        << " cacheAgglomeration:" << cacheAgglomeration_
        << " nPreSweeps:" << nPreSweeps_
        << " preSweepsLevelMultiplier:" << preSweepsLevelMultiplier_
        << " maxPreSweeps:" << maxPreSweeps_
        << " nPostSweeps:" << nPostSweeps_
        << " postSweepsLevelMultiplier:" << postSweepsLevelMultiplier_
        << " maxPostSweeps:" << maxPostSweeps_
        << " nFinestSweeps:" << nFinestSweeps_
        << " interpolateCorrection:" << interpolateCorrection_
        << " scaleCorrection:" << scaleCorrection_
        << " directSolveCoarsest:" << directSolveCoarsest_
        << " processorAgglomerate:" << processorAgglomerate_
        << endl;
}


const Foam::lduMatrix& Foam::GAMGSolver::matrixLevel(const label i) const
{
    if (i == 0)
    {
        return matrix_;
    }
    else
    {
        return matrixLevels_[i - 1];
    }
}


const Foam::lduInterfaceFieldPtrsList& Foam::GAMGSolver::interfaceLevel
(
    const label i
) const
{
    if (i == 0)
    {
        return interfaces_;
    }
    else
    {
        return interfaceLevels_[i - 1];
    }
}


const Foam::FieldField<Foam::Field, Foam::scalar>&
Foam::GAMGSolver::interfaceBouCoeffsLevel
(
    const label i
) const
{
    if (i == 0)
    {
        return interfaceBouCoeffs_;
    }
    else
    {
        return interfaceLevelsBouCoeffs_[i - 1];
    }
}


const Foam::FieldField<Foam::Field, Foam::scalar>&
Foam::GAMGSolver::interfaceIntCoeffsLevel
(
    const label i
) const
{
    if (i == 0)
    {
        return interfaceIntCoeffs_;
    }
    else
    {
        return interfaceLevelsIntCoeffs_[i - 1];
    }
}


// ************************************************************************* //
