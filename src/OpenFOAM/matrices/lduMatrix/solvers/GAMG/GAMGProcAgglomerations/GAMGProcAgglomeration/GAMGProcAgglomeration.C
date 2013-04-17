/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "GAMGProcAgglomeration.H"
#include "GAMGAgglomeration.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GAMGProcAgglomeration, 0);
    defineRunTimeSelectionTable(GAMGProcAgglomeration, GAMGAgglomeration);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::GAMGProcAgglomeration::printStats
(
    Ostream& os,
    GAMGAgglomeration& agglom
) const
{
    for (label levelI = 0; levelI <= agglom.size(); levelI++)
    {
        if (agglom.hasMeshLevel(levelI))
        {
            const lduMesh& fineMesh = agglom.meshLevel(levelI);
            const lduInterfacePtrsList& interfaces =
                agglom.interfaceLevel(levelI);

            os  << "Level " << levelI << " fine mesh:"<< nl
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
                    os  << "        " << i
                        << "\tsize:" << interfaces[i].faceCells().size()
                        << endl;
                }
            }

            os  << fineMesh.info() << endl;

            os  << endl;
        }
        else
        {
            os  << "Level " << levelI << " has no fine mesh:" << nl
                << endl;
        }

        if
        (
            levelI < agglom.restrictAddressing_.size()
         && agglom.restrictAddressing_.set(levelI)
        )
        {
            const labelList& cellRestrict =
                agglom.restrictAddressing(levelI);
            const labelList& faceRestrict =
                agglom.faceRestrictAddressing(levelI);

            os  << "Level " << levelI << " agglomeration:" << nl
                << "    nCoarseCells:" << agglom.nCells(levelI) << nl
                << "    nCoarseFaces:" << agglom.nFaces(levelI) << nl
                << "    cellRestriction:"
                << " size:" << cellRestrict.size()
                << " max:" << max(cellRestrict)
                << nl
                << "    faceRestriction:"
                << " size:" << faceRestrict.size()
                << " max:" << max(faceRestrict)
                << nl;

            const labelListList& patchFaceRestrict =
                agglom.patchFaceRestrictAddressing(levelI);
            forAll(patchFaceRestrict, i)
            {
                if (patchFaceRestrict[i].size())
                {
                    const labelList& faceRestrict =
                        patchFaceRestrict[i];
                    os  << "        " << i
                        << " size:" << faceRestrict.size()
                        << " max:" << max(faceRestrict)
                        << nl;
                }
            }
        }
        if
        (
            levelI < agglom.procCellOffsets_.size()
         && agglom.procCellOffsets_.set(levelI)
        )
        {
            os  << "    procCellOffsets:" << agglom.procCellOffsets_[levelI]
                << nl
                << "    procAgglomMap:" << agglom.procAgglomMap_[levelI]
                << nl
                << "    procIDs:" << agglom.agglomProcIDs_[levelI]
                << nl
                << "    comm:" << agglom.procCommunicator_[levelI]
                << endl;
        }

        os  << endl;
    }
    os  << endl;
}


bool Foam::GAMGProcAgglomeration::agglomerate
(
    const label fineLevelIndex,
    const labelList& procAgglomMap,
    const labelList& masterProcs,
    const List<int>& agglomProcIDs,
    const label procAgglomComm
)
{
    const lduMesh& levelMesh = agglom_.meshLevels_[fineLevelIndex];
    label levelComm = levelMesh.comm();

    if (Pstream::myProcNo(levelComm) != -1)
    {
        // Collect meshes and restrictAddressing onto master
        // Overwrites the fine mesh (meshLevels_[index-1]) and addressing
        // from fine mesh to coarse mesh (restrictAddressing_[index]).
        agglom_.procAgglomerateLduAddressing
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
            levelI < agglom_.meshLevels_.size();
            levelI++
        )
        {
            agglom_.procAgglomerateRestrictAddressing
            (
                levelComm,
                agglomProcIDs,
                levelI
            );
        }

        if (Pstream::myProcNo(levelComm) == agglomProcIDs[0])
        {
            // On master. Recreate coarse meshes from restrict addressing
            for
            (
                label levelI = fineLevelIndex;
                levelI < agglom_.meshLevels_.size();
                levelI++
            )
            {
                agglom_.agglomerateLduAddressing(levelI);
            }
        }
        else
        {
            // Agglomerated away. Clear mesh storage.
            for
            (
                label levelI = fineLevelIndex+1;
                levelI <= agglom_.size();
                levelI++
            )
            {
                agglom_.clearLevel(levelI);
            }
        }
    }

    // Should check!
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GAMGProcAgglomeration::GAMGProcAgglomeration
(
    GAMGAgglomeration& agglom,
    const dictionary& controlDict
)
:
    agglom_(agglom)
{}


Foam::autoPtr<Foam::GAMGProcAgglomeration> Foam::GAMGProcAgglomeration::New
(
    const word& type,
    GAMGAgglomeration& agglom,
    const dictionary& controlDict
)
{
    if (debug)
    {
        Info<< "GAMGProcAgglomeration::New(const word&, GAMGAgglomeration&"
               ", const dictionary&) : "
               "constructing GAMGProcAgglomeration"
            << endl;
    }

    GAMGAgglomerationConstructorTable::iterator cstrIter =
        GAMGAgglomerationConstructorTablePtr_->find(type);

    if (cstrIter == GAMGAgglomerationConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "GAMGProcAgglomeration::New(const word&, GAMGAgglomeration&"
            ", const dictionary&) "
        )   << "Unknown GAMGProcAgglomeration type "
            << type << " for GAMGAgglomeration " << agglom.type() << nl << nl
            << "Valid GAMGProcAgglomeration types are :" << endl
            << GAMGAgglomerationConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<GAMGProcAgglomeration>(cstrIter()(agglom, controlDict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GAMGProcAgglomeration::~GAMGProcAgglomeration()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
