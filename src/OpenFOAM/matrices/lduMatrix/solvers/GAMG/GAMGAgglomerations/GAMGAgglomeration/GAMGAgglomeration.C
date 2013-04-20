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
#include "GAMGProcAgglomeration.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GAMGAgglomeration, 0);
    defineRunTimeSelectionTable(GAMGAgglomeration, lduMesh);
    defineRunTimeSelectionTable(GAMGAgglomeration, lduMatrix);
    defineRunTimeSelectionTable(GAMGAgglomeration, geometry);
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
    if (processorAgglomerate())
    {
        procAgglomMap_.setSize(nCreatedLevels);
        agglomProcIDs_.setSize(nCreatedLevels);
        procCellOffsets_.setSize(nCreatedLevels);
        procFaceMap_.setSize(nCreatedLevels);
        procBoundaryMap_.setSize(nCreatedLevels);
        procBoundaryFaceMap_.setSize(nCreatedLevels);

        procAgglomeratorPtr_().agglomerate();
    }

    if (debug)
    {
        for (label levelI = 0; levelI <= size(); levelI++)
        {
            if (hasMeshLevel(levelI))
            {
                const lduMesh& fineMesh = meshLevel(levelI);
                Pout<< "Level " << levelI << " fine mesh:"<< nl;
                Pout<< fineMesh.info() << endl;
            }
            else
            {
                Pout<< "Level " << levelI << " has no fine mesh:" << nl
                    << endl;
            }
        }
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
    procAgglomeratorPtr_
    (
        (
            (UPstream::nProcs(mesh.comm()) > 1)
         && controlDict.found("processorAgglomerator")
        )
      ? GAMGProcAgglomeration::New
        (
            controlDict.lookup("processorAgglomerator"),
            *this,
            controlDict
        )
      : autoPtr<GAMGProcAgglomeration>(NULL)
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
    if (processorAgglomerate())
    {
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


Foam::autoPtr<Foam::GAMGAgglomeration> Foam::GAMGAgglomeration::New
(
    const lduMesh& mesh,
    const scalarField& cellVolumes,
    const vectorField& faceAreas,
    const dictionary& controlDict
)
{
    const word agglomeratorType(controlDict.lookup("agglomerator"));

    const_cast<Time&>(mesh.thisDb().time()).libs().open
    (
        controlDict,
        "geometricGAMGAgglomerationLibs",
        geometryConstructorTablePtr_
    );

    geometryConstructorTable::iterator cstrIter =
        geometryConstructorTablePtr_->find(agglomeratorType);

    if (cstrIter == geometryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "GAMGAgglomeration::New"
            "(const lduMesh& mesh, const dictionary& controlDict)"
        )   << "Unknown GAMGAgglomeration type "
            << agglomeratorType << ".\n"
            << "Valid geometric GAMGAgglomeration types are :"
            << geometryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<GAMGAgglomeration>
    (
        cstrIter()
        (
            mesh,
            cellVolumes,
            faceAreas,
            controlDict
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GAMGAgglomeration::~GAMGAgglomeration()
{}


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
