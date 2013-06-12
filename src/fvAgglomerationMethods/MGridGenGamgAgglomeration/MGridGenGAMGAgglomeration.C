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

#include "MGridGenGAMGAgglomeration.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MGridGenGAMGAgglomeration, 0);

    addToRunTimeSelectionTable
    (
        GAMGAgglomeration,
        MGridGenGAMGAgglomeration,
        lduMesh
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MGridGenGAMGAgglomeration::MGridGenGAMGAgglomeration
(
    const lduMesh& mesh,
    const dictionary& controlDict
)
:
    GAMGAgglomeration(mesh, controlDict),
    fvMesh_(refCast<const fvMesh>(mesh))
{
    // Min, max size of agglomerated cells
    label minSize(readLabel(controlDict.lookup("minSize")));
    label maxSize(readLabel(controlDict.lookup("maxSize")));

    // Start geometric agglomeration from the cell volumes and areas of the mesh
    scalarField* VPtr = const_cast<scalarField*>(&fvMesh_.cellVolumes());

    vectorField magFaceAreas(vector::one*mag(fvMesh_.faceAreas()));
    SubField<vector> Sf(magFaceAreas, fvMesh_.nInternalFaces());
    //SubField<vector> Sf(fvMesh_.faceAreas(), fvMesh_.nInternalFaces());

    vectorField* SfPtr = const_cast<vectorField*>
    (
        &Sf.operator const vectorField&()
    );

    // Create the boundary area cell field
    scalarField* SbPtr(new scalarField(fvMesh_.nCells(), 0));

    {
        scalarField& Sb = *SbPtr;

        const labelList& own = fvMesh_.faceOwner();
        const vectorField& Sf = fvMesh_.faceAreas();

        forAll(Sf, facei)
        {
            if (!fvMesh_.isInternalFace(facei))
            {
                Sb[own[facei]] += mag(Sf[facei]);
            }
        }
    }


    // Agglomerate until the required number of cells in the coarsest level
    // is reached

    label nCreatedLevels = 0;

    while (nCreatedLevels < maxLevels_ - 1)
    {
        label nCoarseCells = -1;

        tmp<labelField> finalAgglomPtr = agglomerate
        (
            nCoarseCells,
            minSize,
            maxSize,
            meshLevel(nCreatedLevels).lduAddr(),
            *VPtr,
            *SfPtr,
            *SbPtr
        );

        if (continueAgglomerating(nCoarseCells))
        {
            nCells_[nCreatedLevels] = nCoarseCells;
            restrictAddressing_.set(nCreatedLevels, finalAgglomPtr);
        }
        else
        {
            break;
        }

        agglomerateLduAddressing(nCreatedLevels);

        // Agglomerate the cell volumes field for the next level
        {
            scalarField* aggVPtr
            (
                new scalarField(meshLevels_[nCreatedLevels].size())
            );

            // Restrict but no parallel agglomeration (not supported)
            restrictField(*aggVPtr, *VPtr, nCreatedLevels, false);

            if (nCreatedLevels)
            {
                delete VPtr;
            }

            VPtr = aggVPtr;
        }

        // Agglomerate the face areas field for the next level
        {
            vectorField* aggSfPtr
            (
                new vectorField
                (
                    meshLevels_[nCreatedLevels].upperAddr().size(),
                    vector::zero
                )
            );

            restrictFaceField(*aggSfPtr, *SfPtr, nCreatedLevels);

            if (nCreatedLevels)
            {
                delete SfPtr;
            }

            SfPtr = aggSfPtr;
        }

        // Agglomerate the cell boundary areas field for the next level
        {
            scalarField* aggSbPtr
            (
                new scalarField(meshLevels_[nCreatedLevels].size())
            );

            // Restrict but no parallel agglomeration (not supported)
            restrictField(*aggSbPtr, *SbPtr, nCreatedLevels, false);

            delete SbPtr;
            SbPtr = aggSbPtr;
        }

        nCreatedLevels++;
    }

    // Shrink the storage of the levels to those created
    compactLevels(nCreatedLevels);

    // Delete temporary geometry storage
    if (nCreatedLevels)
    {
        delete VPtr;
        delete SfPtr;
    }
    delete SbPtr;
}


// ************************************************************************* //
