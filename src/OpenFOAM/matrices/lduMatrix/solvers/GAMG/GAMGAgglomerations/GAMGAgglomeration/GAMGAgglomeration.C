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

#include "GAMGAgglomeration.H"
#include "lduMesh.H"
#include "lduMatrix.H"
#include "Time.H"
#include "dlLibraryTable.H"
#include "GAMGInterface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GAMGAgglomeration, 0);
    defineRunTimeSelectionTable(GAMGAgglomeration, lduMesh);
    defineRunTimeSelectionTable(GAMGAgglomeration, lduMatrix);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::GAMGAgglomeration::compactLevels(const label nCreatedLevels)
{
    nCells_.setSize(nCreatedLevels);
    restrictAddressing_.setSize(nCreatedLevels);
    nFaces_.setSize(nCreatedLevels);
    faceRestrictAddressing_.setSize(nCreatedLevels);
    faceFlipMap_.setSize(nCreatedLevels);
    nPatchFaces_.setSize(nCreatedLevels);
    patchFaceRestrictAddressing_.setSize(nCreatedLevels);
    meshLevels_.setSize(nCreatedLevels);
    primitiveInterfaces_.setSize(nCreatedLevels + 1);
    interfaceLevels_.setSize(nCreatedLevels + 1);

    // Have procCommunicator_ always, even if not procAgglomerating
    procCommunicator_.setSize(nCreatedLevels + 1);
    if (processorAgglomerate_)
    {
        Pout<< "GAMGAgglomeration : compacting from " << procAgglomMap_.size()
            << " to " << nCreatedLevels << " levels" << endl;
        procAgglomMap_.setSize(nCreatedLevels);
        agglomProcIDs_.setSize(nCreatedLevels);
        procCellOffsets_.setSize(nCreatedLevels);
        procFaceMap_.setSize(nCreatedLevels);
        procBoundaryMap_.setSize(nCreatedLevels);
        procBoundaryFaceMap_.setSize(nCreatedLevels);


        if (debug)
        {
            Pout<< nl << "Starting mesh overview" << endl;
            for (label levelI = 0; levelI <= size(); levelI++)
            {
                if (hasMeshLevel(levelI))
                {
                    const lduMesh& fineMesh = meshLevel(levelI);
                    const lduInterfacePtrsList& interfaces =
                        interfaceLevel(levelI);

                    Pout<< "Level " << levelI << " fine mesh:"<< nl
                        << "    nCells:"
                        << fineMesh.lduAddr().size() << nl
                        << "    nFaces:"
                        << fineMesh.lduAddr().lowerAddr().size() << nl
                        << "    nInterfaces:" << interfaces.size()
                        << endl;

                    forAll(interfaces, i)
                    {
                        if (interfaces.set(i))
                        {
                            Pout<< "        " << i
                                << "\tsize:" << interfaces[i].faceCells().size()
                                << endl;
                        }
                    }
                    Pout<< endl;
                }
                else
                {
                    Pout<< "Level " << levelI << " has no fine mesh:" << nl
                        << endl;
                }

                if
                (
                    levelI < restrictAddressing_.size()
                 && restrictAddressing_.set(levelI)
                )
                {
                    const labelList& cellRestrict = restrictAddressing(levelI);
                    const labelList& faceRestrict =
                        faceRestrictAddressing(levelI);

                    Pout<< "Level " << levelI << " agglomeration:" << nl
                        << "    nCoarseCells:" << nCells(levelI) << nl
                        << "    nCoarseFaces:" << nFaces(levelI) << nl
                        << "    cellRestriction:"
                        << " size:" << cellRestrict.size()
                        << " max:" << max(cellRestrict)
                        << nl
                        << "    faceRestriction:"
                        << " size:" << faceRestrict.size()
                        << " max:" << max(faceRestrict)
                        << nl;


                    const labelListList& patchFaceRestrict =
                        patchFaceRestrictAddressing(levelI);
                    forAll(patchFaceRestrict, i)
                    {
                        if (patchFaceRestrict[i].size())
                        {
                            const labelList& faceRestrict =
                                patchFaceRestrict[i];
                            Pout<< "        " << i
                                << " size:" << faceRestrict.size()
                                << " max:" << max(faceRestrict)
                                << nl;
                        }
                    }
                    Pout<< endl;
                }
                Pout<< endl;
            }
            Pout<< endl;
        }




//XXXXXX
        // As a test: agglomerate finelevel 2, coarselevel 3
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        label fineLevelIndex = 2;

        // Get the coarse mesh
        const lduMesh& levelMesh = meshLevels_[fineLevelIndex];
        label levelComm = levelMesh.comm();

        // Processor restriction map: per processor the coarse processor
        labelList procAgglomMap(UPstream::nProcs(levelComm));
        {
            label half = (procAgglomMap.size()+1)/2;
            for (label i = 0; i < half; i++)
            {
                procAgglomMap[i] = 0;
            }
            for (label i = half; i < procAgglomMap.size(); i++)
            {
                procAgglomMap[i] = 1;
            }
        }

        // Master processor
        labelList masterProcs;
        // Local processors that agglomerate. agglomProcIDs[0] is in
        // masterProc.
        List<int> agglomProcIDs;
        calculateRegionMaster
        (
            levelComm,
            procAgglomMap,
            masterProcs,
            agglomProcIDs
        );

        Pout<< "procAgglomMap:" << procAgglomMap << endl;
        Pout<< "agglomProcIDs:" << agglomProcIDs << endl;
        Pout<< "masterProcs:" << masterProcs << endl;


        // Allocate a communicator for the processor-agglomerated matrix
        label procAgglomComm = UPstream::allocateCommunicator
        (
            levelComm,
            masterProcs
        );


        // Above should really be determined by a procesor-agglomeration
        // engine. Use procesor agglomeration maps to do the actual
        // collecting.

        if (Pstream::myProcNo(levelComm) != -1)
        {
            // Collect meshes and restrictAddressing onto master
            // Overwrites the fine mesh (meshLevels_[index-1]) and addressing
            // from fine mesh to coarse mesh (restrictAddressing_[index]).
            procAgglomerateLduAddressing
            (
                levelComm,
                procAgglomMap,
                agglomProcIDs,
                procAgglomComm,

                fineLevelIndex               //fine level index
            );

            // Combine restrict addressing only onto master
            for
            (
                label levelI = fineLevelIndex+1;
                levelI < meshLevels_.size();
                levelI++
            )
            {
                Pout<< "Starting procAgglomerateRestrictAddressing level:"
                    << levelI << endl;

                procAgglomerateRestrictAddressing
                (
                    levelComm,
                    agglomProcIDs,
                    levelI
                );
                Pout<< "Finished procAgglomerateRestrictAddressing level:"
                    << levelI << endl;
            }

            if (Pstream::myProcNo(levelComm) == agglomProcIDs[0])
            {
                // On master. Recreate coarse meshes from restrict addressing
                for
                (
                    label levelI = fineLevelIndex;
                    levelI < meshLevels_.size();
                    levelI++
                )
                {
                    Pout<< "Starting agglomerateLduAddressing level:" << levelI
                        << endl;
                    agglomerateLduAddressing(levelI);
                    Pout<< "Finished agglomerateLduAddressing level:" << levelI
                        << endl;
                }
            }
            else
            {
                // Agglomerated away. Clear mesh storage.
                for
                (
                    label levelI = fineLevelIndex+1;
                    levelI <= size();
                    levelI++
                )
                {
                    clearLevel(levelI);
                }
            }
        }
    }

    // Print a bit
    if (debug)
    {
        Pout<< nl << "Mesh overview" << endl;
        for (label levelI = 0; levelI <= size(); levelI++)
        {
            if (hasMeshLevel(levelI))
            {
                const lduMesh& fineMesh = meshLevel(levelI);
                const lduInterfacePtrsList& interfaces =
                    interfaceLevel(levelI);

                Pout<< "Level " << levelI << " fine mesh:"<< nl
                    << "    nCells:"
                    << fineMesh.lduAddr().size() << nl
                    << "    nFaces:"
                    << fineMesh.lduAddr().lowerAddr().size() << nl
                    << "    nInterfaces:" << interfaces.size()
                    << endl;

                forAll(interfaces, i)
                {
                    if (interfaces.set(i))
                    {
                        Pout<< "        " << i
                            << "\tsize:" << interfaces[i].faceCells().size()
                            << endl;
                    }
                }

                Pout<< fineMesh.info() << endl;

                Pout<< endl;
            }
            else
            {
                Pout<< "Level " << levelI << " has no fine mesh:" << nl
                    << endl;
            }

            if
            (
                levelI < restrictAddressing_.size()
             && restrictAddressing_.set(levelI)
            )
            {
                const labelList& cellRestrict = restrictAddressing(levelI);
                const labelList& faceRestrict =
                    faceRestrictAddressing(levelI);

                Pout<< "Level " << levelI << " agglomeration:" << nl
                    << "    nCoarseCells:" << nCells(levelI) << nl
                    << "    nCoarseFaces:" << nFaces(levelI) << nl
                    << "    cellRestriction:"
                    << " size:" << cellRestrict.size()
                    << " max:" << max(cellRestrict)
                    << nl
                    << "    faceRestriction:"
                    << " size:" << faceRestrict.size()
                    << " max:" << max(faceRestrict)
                    << nl;

                const labelListList& patchFaceRestrict =
                    patchFaceRestrictAddressing(levelI);
                forAll(patchFaceRestrict, i)
                {
                    if (patchFaceRestrict[i].size())
                    {
                        const labelList& faceRestrict =
                            patchFaceRestrict[i];
                        Pout<< "        " << i
                            << " size:" << faceRestrict.size()
                            << " max:" << max(faceRestrict)
                            << nl;
                    }
                }
            }
            if
            (
                levelI < procCellOffsets_.size()
             && procCellOffsets_.set(levelI)
            )
            {
                Pout<< "    procCellOffsets:" << procCellOffsets_[levelI]
                    << nl
                    << "    procAgglomMap:" << procAgglomMap_[levelI]
                    << nl
                    << "    procIDs:" << agglomProcIDs_[levelI]
                    << nl
                    << "    comm:" << procCommunicator_[levelI]
                    << endl;
            }

            Pout<< endl;
        }
        Pout<< endl;
    }
}


bool Foam::GAMGAgglomeration::continueAgglomerating
(
    const label nCoarseCells
) const
{
    // Check the need for further agglomeration on all processors
    bool contAgg = nCoarseCells >= nCellsInCoarsestLevel_;
    mesh().reduce(contAgg, andOp<bool>());
    return contAgg;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GAMGAgglomeration::GAMGAgglomeration
(
    const lduMesh& mesh,
    const dictionary& controlDict
)
:
    MeshObject<lduMesh, GAMGAgglomeration>(mesh),

    maxLevels_(50),

    nCellsInCoarsestLevel_
    (
        readLabel(controlDict.lookup("nCellsInCoarsestLevel"))
    ),
    processorAgglomerate_
    (
        controlDict.lookupOrDefault<bool>("processorAgglomerate", false)
    ),

    nCells_(maxLevels_),
    restrictAddressing_(maxLevels_),
    nFaces_(maxLevels_),
    faceRestrictAddressing_(maxLevels_),
    faceFlipMap_(maxLevels_),
    nPatchFaces_(maxLevels_),
    patchFaceRestrictAddressing_(maxLevels_),

    meshLevels_(maxLevels_),
    primitiveInterfaces_(maxLevels_ + 1),
    interfaceLevels_(maxLevels_ + 1)
{
    procCommunicator_.setSize(maxLevels_ + 1, -1);
    if (processorAgglomerate_)
    {
        Pout<< "GAMGAgglomeration : sizing to " << maxLevels_
            << " levels" << endl;
        procAgglomMap_.setSize(maxLevels_);
        agglomProcIDs_.setSize(maxLevels_);
        procCellOffsets_.setSize(maxLevels_);
        procFaceMap_.setSize(maxLevels_);
        procBoundaryMap_.setSize(maxLevels_);
        procBoundaryFaceMap_.setSize(maxLevels_);
    }
}


const Foam::GAMGAgglomeration& Foam::GAMGAgglomeration::New
(
    const lduMesh& mesh,
    const dictionary& controlDict
)
{
    if
    (
        !mesh.thisDb().foundObject<GAMGAgglomeration>
        (
            GAMGAgglomeration::typeName
        )
    )
    {
        const word agglomeratorType(controlDict.lookup("agglomerator"));

        const_cast<Time&>(mesh.thisDb().time()).libs().open
        (
            controlDict,
            "geometricGAMGAgglomerationLibs",
            lduMeshConstructorTablePtr_
        );

        lduMeshConstructorTable::iterator cstrIter =
            lduMeshConstructorTablePtr_->find(agglomeratorType);

        if (cstrIter == lduMeshConstructorTablePtr_->end())
        {
            FatalErrorIn
            (
                "GAMGAgglomeration::New"
                "(const lduMesh& mesh, const dictionary& controlDict)"
            )   << "Unknown GAMGAgglomeration type "
                << agglomeratorType << ".\n"
                << "Valid algebraic GAMGAgglomeration types are :"
                << lduMatrixConstructorTablePtr_->sortedToc() << endl
                << "Valid algebraic GAMGAgglomeration types are :"
                << lduMeshConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return store(cstrIter()(mesh, controlDict).ptr());
    }
    else
    {
        return mesh.thisDb().lookupObject<GAMGAgglomeration>
        (
            GAMGAgglomeration::typeName
        );
    }
}


const Foam::GAMGAgglomeration& Foam::GAMGAgglomeration::New
(
    const lduMatrix& matrix,
    const dictionary& controlDict
)
{
    const lduMesh& mesh = matrix.mesh();

    if
    (
        !mesh.thisDb().foundObject<GAMGAgglomeration>
        (
            GAMGAgglomeration::typeName
        )
    )
    {
        const word agglomeratorType(controlDict.lookup("agglomerator"));

        const_cast<Time&>(mesh.thisDb().time()).libs().open
        (
            controlDict,
            "algebraicGAMGAgglomerationLibs",
            lduMatrixConstructorTablePtr_
        );

        if
        (
            !lduMatrixConstructorTablePtr_
         || !lduMatrixConstructorTablePtr_->found(agglomeratorType)
        )
        {
            return New(mesh, controlDict);
        }
        else
        {
            lduMatrixConstructorTable::iterator cstrIter =
                lduMatrixConstructorTablePtr_->find(agglomeratorType);

            return store(cstrIter()(matrix, controlDict).ptr());
        }
    }
    else
    {
        return mesh.thisDb().lookupObject<GAMGAgglomeration>
        (
            GAMGAgglomeration::typeName
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GAMGAgglomeration::~GAMGAgglomeration()
{
    forAllReverse(procCommunicator_, i)
    {
        if (procCommunicator_[i] != -1)
        {
            UPstream::freeCommunicator(procCommunicator_[i]);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::lduMesh& Foam::GAMGAgglomeration::meshLevel
(
    const label i
) const
{
    if (i == 0)
    {
        return mesh_;
    }
    else
    {
        return meshLevels_[i - 1];
    }
}


bool Foam::GAMGAgglomeration::hasMeshLevel(const label i) const
{
    if (i == 0)
    {
        return true;
    }
    else
    {
        return meshLevels_.set(i - 1);
    }
}


const Foam::lduInterfacePtrsList& Foam::GAMGAgglomeration::interfaceLevel
(
    const label i
) const
{
    return interfaceLevels_[i];
}


void Foam::GAMGAgglomeration::clearLevel(const label i)
{
    if (hasMeshLevel(i))
    {
        Pout<< "Clearing out level " << i << endl;

        meshLevels_.set(i - 1, NULL);

        if (i < nCells_.size())
        {
            nCells_[i] = -555;
            restrictAddressing_.set(i, NULL);
            nFaces_[i] = -666;
            faceRestrictAddressing_.set(i, NULL);
            faceFlipMap_.set(i, NULL);
            nPatchFaces_.set(i, NULL);
            patchFaceRestrictAddressing_.set(i, NULL);
            primitiveInterfaces_.set(i, NULL);
            interfaceLevels_.set(i, NULL);
        }
    }
}


const Foam::labelList& Foam::GAMGAgglomeration::procAgglomMap
(
    const label leveli
) const
{
    return procAgglomMap_[leveli];
}


const Foam::labelList& Foam::GAMGAgglomeration::agglomProcIDs
(
    const label leveli
) const
{
    return agglomProcIDs_[leveli];
}


bool Foam::GAMGAgglomeration::hasProcMesh(const label leveli) const
{
    return procCommunicator_[leveli] != -1;
}


Foam::label Foam::GAMGAgglomeration::procCommunicator(const label leveli) const
{
    return procCommunicator_[leveli];
}


const Foam::labelList& Foam::GAMGAgglomeration::cellOffsets
(
    const label leveli
) const
{
    return procCellOffsets_[leveli];
}


const Foam::labelListList& Foam::GAMGAgglomeration::faceMap
(
    const label leveli
) const
{
    return procFaceMap_[leveli];
}


const Foam::labelListList& Foam::GAMGAgglomeration::boundaryMap
(
    const label leveli
) const
{
    return procBoundaryMap_[leveli];
}


const Foam::labelListListList& Foam::GAMGAgglomeration::boundaryFaceMap
(
    const label leveli
) const
{
    return procBoundaryFaceMap_[leveli];
}


// ************************************************************************* //
