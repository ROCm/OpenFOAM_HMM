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
#include "GAMGInterface.H"

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
    agglomeration_(GAMGAgglomeration::New(matrix_, controlDict_)),

    matrixLevels_(agglomeration_.size()),
    interfaceLevels_(agglomeration_.size()),
    primitiveInterfaceLevels_(agglomeration_.size()),
    interfaceLevelsBouCoeffs_(agglomeration_.size()),
    interfaceLevelsIntCoeffs_(agglomeration_.size())
{
    readControls();

    //if (agglomeration_.processorAgglomerate())
    //{
    //    procMatrixLevels_.setSize(agglomeration_.size());
    //    procInterfaceLevels_.setSize(agglomeration_.size());
    //    procPrimitiveInterfaces_.setSize(agglomeration_.size());
    //    procInterfaceLevelsBouCoeffs_.setSize(agglomeration_.size());
    //    procInterfaceLevelsIntCoeffs_.setSize(agglomeration_.size());
    //}

    if (agglomeration_.processorAgglomerate())
    {
        forAll(agglomeration_, fineLevelIndex)
        {
            if (agglomeration_.hasMeshLevel(fineLevelIndex))
            {
                Pout<< "Level:" << fineLevelIndex
                    << " agglomerating matrix." << endl;

                if
                (
                    //   !agglomeration_.hasMeshLevel(fineLevelIndex+1)
                    //||  (
                    //        agglomeration_.nCells(fineLevelIndex)
                    //     != agglomeration_.meshLevel(fineLevelIndex+1).
                    //            lduAddr().size()
                    //    )

                    (fineLevelIndex+1) < agglomeration_.procBoundaryMap_.size()
                 && agglomeration_.procBoundaryMap_.set(fineLevelIndex+1)
                )
                {
                    Pout<< "Level:" << fineLevelIndex
                        << " agglomerating onto dummy coarse mesh." << endl;

                    // Construct matrix without referencing the coarse mesh so
                    // construct a dummy mesh instead. This will get overwritten
                    // by the call to procAgglomerateMatrix so is only to get
                    // it through agglomerateMatrix


                    const lduInterfacePtrsList& fineMeshInterfaces =
                        agglomeration_.interfaceLevel(fineLevelIndex);

                    PtrList<GAMGInterface> dummyPrimMeshInterfaces
                    (
                        fineMeshInterfaces.size()
                    );
                    lduInterfacePtrsList dummyMeshInterfaces
                    (
                        dummyPrimMeshInterfaces.size()
                    );
                    forAll(fineMeshInterfaces, intI)
                    {
                        if (fineMeshInterfaces.set(intI))
                        {
                            OStringStream os;
                            refCast<const GAMGInterface>
                            (
                                fineMeshInterfaces[intI]
                            ).write(os);
                            IStringStream is(os.str());

                            dummyPrimMeshInterfaces.set
                            (
                                intI,
                                GAMGInterface::New
                                (
                                    fineMeshInterfaces[intI].type(),
                                    intI,
                                    dummyMeshInterfaces,
                                    is
                                )
                            );
                        }
                    }

                    forAll(dummyPrimMeshInterfaces, intI)
                    {
                        if (dummyPrimMeshInterfaces.set(intI))
                        {
                            dummyMeshInterfaces.set
                            (
                                intI,
                                &dummyPrimMeshInterfaces[intI]
                            );
                        }
                    }

                    // So:
                    // - pass in incorrect mesh (= fine mesh instead of coarse)
                    // - pass in dummy interfaces
                    agglomerateMatrix
                    (
                        fineLevelIndex,
                        agglomeration_.meshLevel(fineLevelIndex),
                        dummyMeshInterfaces
                    );


                    const labelList& procAgglomMap =
                        agglomeration_.procAgglomMap(fineLevelIndex+1);
                    const List<int>& procIDs =
                        agglomeration_.agglomProcIDs(fineLevelIndex+1);

                    Pout<< "Level:" << fineLevelIndex
                        << " agglomerating onto procIDs:" << procIDs << endl;

                    procAgglomerateMatrix
                    (
                        procAgglomMap,
                        procIDs,
                        fineLevelIndex
                    );

                    Pout<< "Level:" << fineLevelIndex
                        << " DONE agglomerating onto procIDs:" << procIDs
                        << endl;
                }
                else
                {
                    Pout<< "Level:" << fineLevelIndex
                        << " agglomerating onto coarse mesh at level "
                        << fineLevelIndex + 1 << endl;

                    agglomerateMatrix
                    (
                        fineLevelIndex,
                        agglomeration_.meshLevel(fineLevelIndex + 1),
                        agglomeration_.interfaceLevel(fineLevelIndex + 1)
                    );
                    Pout<< "Level:" << fineLevelIndex
                        << " DONE agglomerating onto coarse mesh at level "
                        << fineLevelIndex + 1 << endl;

                }
            }
            else
            {
                // No mesh. Not involved in calculation anymore
            }
        }
    }
    else
    {
        forAll(agglomeration_, fineLevelIndex)
        {
            // Agglomerate on to coarse level mesh
            agglomerateMatrix
            (
                fineLevelIndex,
                agglomeration_.meshLevel(fineLevelIndex + 1),
                agglomeration_.interfaceLevel(fineLevelIndex + 1)
            );
        }
    }


    if (debug)
    {
        for
        (
            label fineLevelIndex = 0;
            fineLevelIndex <= matrixLevels_.size();
            fineLevelIndex++
        )
        {
            if (fineLevelIndex == 0 || matrixLevels_.set(fineLevelIndex-1))
            {
                const lduMatrix& matrix = matrixLevel(fineLevelIndex);
                const lduInterfaceFieldPtrsList& interfaces =
                    interfaceLevel(fineLevelIndex);

                Pout<< "level:" << fineLevelIndex << nl
                    << "    nCells:" << matrix.diag().size() << nl
                    << "    nFaces:" << matrix.lower().size() << nl
                    << "    nInterfaces:" << interfaces.size()
                    << endl;

                forAll(interfaces, i)
                {
                    if (interfaces.set(i))
                    {
                        Pout<< "        " << i
                            << "\ttype:" << interfaces[i].type()
                            << endl;
                    }
                }
            }
            else
            {
                Pout<< "level:" << fineLevelIndex << " : no matrix" << endl;
            }
        }
        Pout<< endl;
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
        else if (agglomeration_.processorAgglomerate())
        {
            //// Pick a level to processor agglomerate
            //label agglomLevel = matrixLevels_.size() - 1;//1;
            //
            //
            //// Get mesh and matrix at this level
            //const lduMatrix& levelMatrix = matrixLevels_[agglomLevel];
            //const lduMesh& levelMesh = levelMatrix.mesh();
            //
            //
            //label levelComm = levelMesh.comm();
            //label oldWarn = UPstream::warnComm;
            //UPstream::warnComm = levelComm;
            //
            //Pout<< "Solve generic on mesh (level=" << agglomLevel
            //    << ") using communicator " << levelComm << endl;
            //
            //// Processor restriction map: per processor the coarse processor
            //labelList procAgglomMap(UPstream::nProcs(levelComm));
            //// Master processor
            //labelList masterProcs;
            //// Local processors that agglomerate. agglomProcIDs[0] is in
            //// masterProc.
            //List<int> agglomProcIDs;
            //
            //{
            //    procAgglomMap[0] = 0;
            //    procAgglomMap[1] = 0;
            //    procAgglomMap[2] = 1;
            //    procAgglomMap[3] = 1;
            //
            //    // Determine the master processors
            //    Map<label> agglomToMaster(procAgglomMap.size());
            //
            //    forAll(procAgglomMap, procI)
            //    {
            //        label coarseI = procAgglomMap[procI];
            //
            //        Map<label>::iterator fnd = agglomToMaster.find(coarseI);
            //        if (fnd == agglomToMaster.end())
            //        {
            //            agglomToMaster.insert(coarseI, procI);
            //        }
            //        else
            //        {
            //            fnd() = max(fnd(), procI);
            //        }
            //    }
            //
            //    masterProcs.setSize(agglomToMaster.size());
            //    forAllConstIter(Map<label>, agglomToMaster, iter)
            //    {
            //        masterProcs[iter.key()] = iter();
            //    }
            //
            //
            //    // Collect all the processors in my agglomeration
            //    label myProcID = Pstream::myProcNo(levelComm);
            //    label myAgglom = procAgglomMap[myProcID];
            //
            //    // Get all processors agglomerating to the same coarse
            //    // processor
            //    agglomProcIDs = findIndices(procAgglomMap, myAgglom);
            //    // Make sure the master is the first element.
            //    label index = findIndex
            //    (
            //        agglomProcIDs,
            //        agglomToMaster[myAgglom]
            //    );
            //    Swap(agglomProcIDs[0], agglomProcIDs[index]);
            //}
            //
            //
            //Pout<< "procAgglomMap:" << procAgglomMap << endl;
            //Pout<< "agglomProcIDs:" << agglomProcIDs << endl;
            //
            //// Allocate a communicator for the processor-agglomerated matrix
            //label procAgglomComm = UPstream::allocateCommunicator
            //(
            //    levelComm,
            //    masterProcs
            //);
            //Pout<< "** Allocated communicator " << procAgglomComm
            //    << " for indices " << masterProcs
            //    << " in processor list " << UPstream::procID(levelComm)
            //    << endl;
            //
            //
            //
            //// Gather matrix and mesh onto agglomProcIDs[0]
            //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //
            //procAgglomerateMatrix
            //(
            //    // Agglomeration information
            //    procAgglomMap,
            //    agglomProcIDs,
            //    procAgglomComm,
            //
            //    agglomLevel,  // level (coarse, not fine level!)
            //
            //    // Resulting matrix
            //    allMatrixPtr_,
            //    allInterfaceBouCoeffs_,
            //    allInterfaceIntCoeffs_,
            //    allPrimitiveInterfaces_,
            //    allInterfaces_
            //);
            //
            //
            //UPstream::warnComm = oldWarn;
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
//    // Clear the the lists of pointers to the interfaces
//    forAll(interfaceLevels_, leveli)
//    {
//        lduInterfaceFieldPtrsList& curLevel = interfaceLevels_[leveli];
//
//        forAll(curLevel, i)
//        {
//            if (curLevel.set(i))
//            {
//                delete curLevel(i);
//            }
//        }
//    }

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
