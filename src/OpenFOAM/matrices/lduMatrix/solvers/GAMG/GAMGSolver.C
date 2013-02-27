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
#include "processorGAMGInterface.H"
#include "processorLduInterfaceField.H"
#include "GAMGInterfaceField.H"
#include "processorGAMGInterfaceField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GAMGSolver, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<GAMGSolver>
        addGAMGSolverMatrixConstructorToTable_;

    lduMatrix::solver::addasymMatrixConstructorToTable<GAMGSolver>
        addGAMGAsymSolverMatrixConstructorToTable_;
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::GAMGSolver::gatherMatrices
(
    const labelList& procIDs,
    const PtrList<lduMesh>& procMeshes,
    const label meshComm,

    const lduMatrix& mat,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,

    PtrList<lduMatrix>& otherMats,
    PtrList<FieldField<Field, scalar> >& otherBouCoeffs,
    PtrList<FieldField<Field, scalar> >& otherIntCoeffs,
    PtrList<lduInterfaceFieldPtrsList>& otherInterfaces
) const
{
    Pout<< "GAMGSolver::gatherMatrices :"
        << " collecting matrices on procs:" << procIDs
        << " using comm:" << meshComm << endl;


    //lduInterfacePtrsList interfaces(mesh.interfaces());

    if (Pstream::myProcNo(meshComm) == procIDs[0])
    {
        // Master.

        Pout<< "GAMGSolver::gatherMatrices :"
            << " master:" << Pstream::myProcNo(meshComm) << endl;

        for (label i=0; i < procMeshes.size(); i++)
        {
            if (meshComm != procMeshes[i].comm())
            {
                FatalErrorIn("GAMGSolver::gatherMatrices()")
                    << "Inconsistent communicator :"
                    << "master processor:" << procIDs[0]
                    << " comm:" << meshComm
                    << " slave processor:" << procIDs[i]
                    << " comm:" << procMeshes[i].comm()
                    << abort(FatalError);
            }
        }


        otherMats.setSize(procIDs.size()-1);
        otherBouCoeffs.setSize(procIDs.size()-1);
        otherIntCoeffs.setSize(procIDs.size()-1);
        otherInterfaces.setSize(procIDs.size()-1);

        for (label i = 1; i < procIDs.size(); i++)
        {
            const lduMesh& procMesh = procMeshes[i];
            const lduInterfacePtrsList procInterfaces = procMesh.interfaces();

            IPstream fromSlave
            (
                Pstream::scheduled,
                procIDs[i],
                0,          // bufSize
                Pstream::msgType(),
                meshComm
            );

            otherMats.set(i-1, new lduMatrix(procMesh, fromSlave));

            // Receive number of/valid interfaces
            boolList validTransforms(fromSlave);
            List<int> validRanks(fromSlave);

            // Size coefficients
            otherBouCoeffs.set
            (
                i-1,
                new FieldField<Field, scalar>(validTransforms.size())
            );
            otherIntCoeffs.set
            (
                i-1,
                new FieldField<Field, scalar>(validTransforms.size())
            );
            otherInterfaces.set
            (
                i-1,
                new lduInterfaceFieldPtrsList(validTransforms.size())
            );

            forAll(validTransforms, intI)
            {
                if (validTransforms[intI])
                {
                    const processorGAMGInterface& interface =
                        refCast<const processorGAMGInterface>
                        (
                            procInterfaces[intI]
                        );


                    otherBouCoeffs[i-1].set(intI, new scalarField(fromSlave));
                    otherIntCoeffs[i-1].set(intI, new scalarField(fromSlave));
                    otherInterfaces[i-1].set
                    (
                        intI,
                        GAMGInterfaceField::New
                        (
                            interface,  //procInterfaces[intI],
                            validTransforms[intI],
                            validRanks[intI]
                        ).ptr()
                    );
                }
            }
        }
    }
    else
    {
        Pout<< "GAMGSolver::gatherMatrices :"
            << " sending from:" << Pstream::myProcNo(meshComm)
            << " to master:" << procIDs[0] << endl;

        // Send to master

        // Count valid interfaces
        boolList validTransforms(interfaceBouCoeffs.size(), false);
        List<int> validRanks(interfaceBouCoeffs.size(), -1);
        forAll(interfaces, intI)
        {
            if (interfaces.set(intI))
            {
                const processorLduInterfaceField& interface =
                    refCast<const processorLduInterfaceField>
                    (
                        interfaces[intI]
                    );

                validTransforms[intI] = interface.doTransform();
                validRanks[intI] = interface.rank();
            }
        }

        Pout<< "GAMGSolver::gatherMatrices :"
            << " sending matrix:" << mat.info() << endl;

        OPstream toMaster
        (
            Pstream::scheduled,
            procIDs[0],
            0,
            Pstream::msgType(),
            meshComm
        );

        toMaster << mat << validTransforms << validRanks;
        forAll(validTransforms, intI)
        {
            if (validTransforms[intI])
            {
                toMaster
                    << interfaceBouCoeffs[intI]
                    << interfaceIntCoeffs[intI];
            }
        }
    }
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
        const label coarsestLevel = matrixLevels_.size() - 1;

        if (directSolveCoarsest_)
        {
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
            const lduMatrix& coarsestMatrix = matrixLevels_[coarsestLevel];
            const lduInterfaceFieldPtrsList& coarsestInterfaces =
                interfaceLevels_[coarsestLevel];
            const FieldField<Field, scalar>& coarsestBouCoeffs =
                interfaceLevelsBouCoeffs_[coarsestLevel];
            const FieldField<Field, scalar>& coarsestIntCoeffs =
                interfaceLevelsIntCoeffs_[coarsestLevel];
            const lduMesh& coarsestMesh = coarsestMatrix.mesh();


            label coarseComm = coarsestMesh.comm();
            label oldWarn = UPstream::warnComm;
            UPstream::warnComm = coarseComm;

            Pout<< "Solve generic on coarsestmesh (level=" << coarsestLevel
                << ") using communicator " << coarseComm << endl;

            // All processors to master
            const List<int>& procIDs = UPstream::procID(coarseComm);

            // Gather all meshes onto procIDs[0]
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // Note: move into GAMGAgglomeration

            Pout<< "procIDs:" << procIDs << endl;

            PtrList<lduMesh> procMeshes;
            lduPrimitiveMesh::gather(coarsestMesh, procIDs, procMeshes);

            forAll(procMeshes, procI)
            {
                Pout<< "** procMesh " << procI << " "
                    << procMeshes[procI].info()
                    << endl;
            }
            Pout<< endl;


            // Gather all matrix coefficients onto procIDs[0]
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            PtrList<lduMatrix> otherMats;
            PtrList<FieldField<Field, scalar> > otherBouCoeffs;
            PtrList<FieldField<Field, scalar> > otherIntCoeffs;
            PtrList<lduInterfaceFieldPtrsList> otherInterfaces;
            gatherMatrices
            (
                procIDs,
                procMeshes,
                coarseComm,

                coarsestMatrix,
                coarsestBouCoeffs,
                coarsestIntCoeffs,
                coarsestInterfaces,

                otherMats,
                otherBouCoeffs,
                otherIntCoeffs,
                otherInterfaces
            );

            Pout<< "Own matrix:" << coarsestMatrix.info() << endl;

            forAll(otherMats, i)
            {
                Pout<< "** otherMats " << i << " "
                    << otherMats[i].info()
                    << endl;
            }
            Pout<< endl;


            // Allocate a communicator for single processor
            label masterComm = UPstream::allocateCommunicator
            (
                coarseComm,
                labelList(1, 0)
            );
            Pout<< "** Allocated communicator " << masterComm
                << " for processor " << procIDs[0]
                << " out of " << procIDs << endl;


            if (Pstream::myProcNo(coarseComm) == procIDs[0])
            {
                // Agglomerate all addressing
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~

                allMeshPtr_.reset
                (
                    new lduPrimitiveMesh
                    (
                        masterComm,
                        procIDs,
                        procMeshes,

                        cellOffsets_,
                        faceMap_,
                        boundaryMap_,
                        boundaryFaceMap_
                    )
                );
                const lduMesh& allMesh = allMeshPtr_();


                // Agglomerate all matrix
                // ~~~~~~~~~~~~~~~~~~~~~~

                allMatrixPtr_.reset(new lduMatrix(allMesh));
                lduMatrix& allMatrix = allMatrixPtr_();

                if (coarsestMatrix.hasDiag())
                {
                    scalarField& allDiag = allMatrix.diag();
                    SubList<scalar>
                    (
                        allDiag,
                        coarsestMatrix.diag().size()
                    ).assign
                    (
                        coarsestMatrix.diag()
                    );
                    forAll(otherMats, i)
                    {
                        SubList<scalar>
                        (
                            allDiag,
                            otherMats[i].diag().size(),
                            cellOffsets_[i+1]
                        ).assign
                        (
                            otherMats[i].diag()
                        );
                    }
                }
                if (coarsestMatrix.hasLower())
                {
                    scalarField& allLower = allMatrix.lower();
                    UIndirectList<scalar>(allLower, faceMap_[0]) =
                        coarsestMatrix.lower();
                    forAll(otherMats, i)
                    {
                        UIndirectList<scalar>(allLower, faceMap_[i+1]) =
                            otherMats[i].lower();
                    }
                }
                if (coarsestMatrix.hasUpper())
                {
                    scalarField& allUpper = allMatrix.upper();
                    UIndirectList<scalar>(allUpper, faceMap_[0]) =
                        coarsestMatrix.upper();
                    forAll(otherMats, i)
                    {
                        UIndirectList<scalar>(allUpper, faceMap_[i+1]) =
                            otherMats[i].upper();
                    }
                }


                // Agglomerate interface fields and coefficients
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                lduInterfacePtrsList allMeshInterfaces = allMesh.interfaces();

                allInterfaceBouCoeffs_.setSize(allMeshInterfaces.size());
                allInterfaceIntCoeffs_.setSize(allMeshInterfaces.size());
                allInterfaces_.setSize(allMeshInterfaces.size());

                forAll(allMeshInterfaces, intI)
                {
                    const lduInterface& patch = allMeshInterfaces[intI];
                    label size = patch.faceCells().size();

                    allInterfaceBouCoeffs_.set(intI, new scalarField(size));
                    allInterfaceIntCoeffs_.set(intI, new scalarField(size));
                }

                labelList nBounFaces(allMeshInterfaces.size());
                forAll(boundaryMap_, procI)
                {
                    const FieldField<Field, scalar>& procBouCoeffs
                    (
                        (procI == 0)
                      ? coarsestBouCoeffs
                      : otherBouCoeffs[procI-1]
                    );
                    const FieldField<Field, scalar>& procIntCoeffs
                    (
                        (procI == 0)
                      ? coarsestIntCoeffs
                      : otherIntCoeffs[procI-1]
                    );
                    const lduInterfaceFieldPtrsList procInterfaces
                    (
                        (procI == 0)
                      ? coarsestInterfaces
                      : otherInterfaces[procI-1]
                    );


                    const labelList& bMap = boundaryMap_[procI];
                    forAll(bMap, procIntI)
                    {
                        label allIntI = bMap[procIntI];

                        if (allIntI != -1)
                        {
                            // So this boundary has been preserved. Copy
                            // data across.

                            if (!allInterfaces_.set(allIntI))
                            {
                                // Construct lduInterfaceField

                                const processorGAMGInterfaceField& procInt =
                                    refCast<const processorGAMGInterfaceField>
                                    (
                                        procInterfaces[procIntI]
                                    );

                                allInterfaces_.set
                                (
                                    allIntI,
                                    GAMGInterfaceField::New
                                    (
                                        refCast<const GAMGInterface>
                                        (
                                            allMeshInterfaces[allIntI]
                                        ),
                                        procInt.doTransform(),
                                        procInt.rank()
                                    ).ptr()
                                );
                            }


                            // Map data from processor to complete mesh

                            scalarField& allBou =
                                allInterfaceBouCoeffs_[allIntI];
                            scalarField& allInt =
                                allInterfaceIntCoeffs_[allIntI];

                            const labelList& map =
                                boundaryFaceMap_[procI][procIntI];
                            const scalarField& procBou =
                                procBouCoeffs[procIntI];
                            const scalarField& procInt =
                                procIntCoeffs[procIntI];

                            forAll(map, i)
                            {
                                label allFaceI = map[i];
                                if (allFaceI < 0)
                                {
                                    FatalErrorIn("GAMGSolver::GAMGSolver()")
                                        << "problem." << abort(FatalError);
                                }
                                allBou[allFaceI] = procBou[i];
                                allInt[allFaceI] = procInt[i];
                            }
                        }
                        else if (procBouCoeffs.set(procIntI))
                        {
                            // Boundary has become internal face

                            // Problem: we don't know whether face is flipped
                            // or not. Either we interrogate the original
                            // mesh (requires procMeshes which we don't want)
                            // or we need to store this information somehow
                            // in the mapping (boundaryMap_, boundaryFaceMap_)


                            //Hack. See above.
                            //const labelList& fCells =
                            //    procInterfaces[procI].faceCells();
                            const labelList& fCells =
                                procMeshes[procI].lduAddr().patchAddr(procIntI);

                            const labelList& map =
                                boundaryFaceMap_[procI][procIntI];

                            const labelList& allU =
                                allMesh.lduAddr().upperAddr();
                            const labelList& allL =
                                allMesh.lduAddr().lowerAddr();
                            const label off = cellOffsets_[procI];

                            const scalarField& procBou =
                                procBouCoeffs[procIntI];
                            const scalarField& procInt =
                                procIntCoeffs[procIntI];

                            forAll(map, i)
                            {
                                if (map[i] >= 0)
                                {
                                    label allFaceI = map[i];

                                    if (coarsestMatrix.hasUpper())
                                    {
                                        allMatrix.upper()[allFaceI] =
                                            procBou[i];
                                    }
                                    if (coarsestMatrix.hasLower())
                                    {
                                        allMatrix.lower()[allFaceI] =
                                            procInt[i];
                                    }
                                }
                                else
                                {
                                    label allFaceI = -map[i]-1;

                                    if (coarsestMatrix.hasUpper())
                                    {
                                        allMatrix.upper()[allFaceI] =
                                            procInt[i];
                                    }
                                    if (coarsestMatrix.hasLower())
                                    {
                                        allMatrix.lower()[allFaceI] =
                                            procBou[i];
                                    }
                                }


                                // Simple check
                                label allFaceI =
                                (
                                    map[i] >= 0
                                  ? map[i]
                                  : -map[i]-1
                                );
                                label allCellI = off + fCells[i];

                                if
                                (
                                    allCellI != allL[allFaceI]
                                 && allCellI != allU[allFaceI]
                                )
                                {
                                    FatalErrorIn
                                    (
                                        "GAMGSolver::GAMGSolver()"
                                    )   << "problem."
                                        << " allFaceI:" << allFaceI
                                        << " local cellI:" << fCells[i]
                                        << " allCellI:" << allCellI
                                        << " allLower:" << allL[allFaceI]
                                        << " allUpper:" << allU[allFaceI]
                                        << abort(FatalError);
                                }
                            }
                        }
                    }
                }
            }

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
