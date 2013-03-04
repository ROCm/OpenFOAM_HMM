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
    meshLevels_.setSize(nCreatedLevels);
    interfaceLevels_.setSize(nCreatedLevels + 1);
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

    nCells_(maxLevels_),
    restrictAddressing_(maxLevels_),
    faceRestrictAddressing_(maxLevels_),

    meshLevels_(maxLevels_),
    interfaceLevels_(maxLevels_ + 1)
{}


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
    // Temporary store the user-defined communicators so we can delete them
    labelList communicators(meshLevels_.size());
    forAll(meshLevels_, leveli)
    {
        communicators[leveli] = meshLevels_[leveli].comm();
    }

    Pout<< "~GAMGAgglomeration() : current communicators:" << communicators
        << endl;

    // Clear the interface storage by hand.
    // It is a list of ptrs not a PtrList for consistency of the interface
    for (label leveli=1; leveli<interfaceLevels_.size(); leveli++)
    {
        lduInterfacePtrsList& curLevel = interfaceLevels_[leveli];

        forAll(curLevel, i)
        {
            if (curLevel.set(i))
            {
                delete curLevel(i);
            }
        }
    }

    forAll(communicators, i)
    {
        UPstream::freeCommunicator(communicators[i]);
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


const Foam::lduInterfacePtrsList& Foam::GAMGAgglomeration::interfaceLevel
(
    const label i
) const
{
    return interfaceLevels_[i];
}


void Foam::GAMGAgglomeration::gatherMeshes
(
    const labelList& procAgglomMap,
    const labelList& procIDs,
    const label allMeshComm
) const
{
    if (allMeshPtr_.valid())
    {
        FatalErrorIn("GAMGAgglomeration::gatherMeshes(..)")
            << "Processor-agglomerated mesh already constructed"
            << exit(FatalError);
    }

    const lduMesh& coarsestMesh = meshLevels_.last();

    label coarseComm = coarsestMesh.comm();

    label oldWarn = UPstream::warnComm;
    UPstream::warnComm = coarseComm;

    Pout<< "GAMGAgglomeration : gathering coarsestmesh (level="
        << meshLevels_.size()-1
        << ") using communicator " << coarseComm << endl;

    lduPrimitiveMesh::gather(coarsestMesh, procIDs, otherMeshes_);

    Pout<< "** Own Mesh " << coarsestMesh.info() << endl;
    forAll(otherMeshes_, i)
    {
        Pout<< "** otherMesh " << i << " "
            << otherMeshes_[i].info()
            << endl;
    }
    Pout<< endl;

    if (Pstream::myProcNo(coarseComm) == procIDs[0])
    {
        // Agglomerate all addressing
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~

        allMeshPtr_.reset
        (
            new lduPrimitiveMesh
            (
                allMeshComm,
                procAgglomMap,

                procIDs,
                coarsestMesh,
                otherMeshes_,

                cellOffsets_,
                faceMap_,
                boundaryMap_,
                boundaryFaceMap_
            )
        );

        Pout<< "** Agglomerated Mesh " << allMeshPtr_().info() << endl;
    }

    UPstream::warnComm = oldWarn;
}


const Foam::lduPrimitiveMesh& Foam::GAMGAgglomeration::allMesh() const
{
    if (!allMeshPtr_.valid())
    {
        FatalErrorIn("GAMGAgglomeration::allMesh() const")
            << "Processor-agglomerated mesh not constructed. Did you call"
            << " gatherMeshes?"
            << exit(FatalError);
    }
    return allMeshPtr_();
}


const Foam::labelList& Foam::GAMGAgglomeration::cellOffsets() const
{
    return cellOffsets_;
}


const Foam::labelListList& Foam::GAMGAgglomeration::faceMap() const
{
    return faceMap_;
}


const Foam::labelListList& Foam::GAMGAgglomeration::boundaryMap() const
{
    return boundaryMap_;
}


const Foam::labelListListList& Foam::GAMGAgglomeration::boundaryFaceMap() const
{
    return boundaryFaceMap_;
}


const Foam::PtrList<Foam::lduMesh>& Foam::GAMGAgglomeration::otherMeshes() const
{
    return otherMeshes_;
}


// ************************************************************************* //
