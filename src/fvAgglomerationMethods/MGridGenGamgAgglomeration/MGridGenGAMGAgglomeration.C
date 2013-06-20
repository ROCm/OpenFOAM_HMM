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
#include "processorPolyPatch.H"
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

    scalarField magFaceAreas(sqrt(3.0)*mag(fvMesh_.faceAreas()));
    SubField<scalar> magSf(magFaceAreas, fvMesh_.nInternalFaces());

    scalarField* magSfPtr = const_cast<scalarField*>
    (
        &magSf.operator const scalarField&()
    );

    // Create the boundary area cell field
    scalarField* magSbPtr(new scalarField(fvMesh_.nCells(), 0));

    {
        scalarField& magSb = *magSbPtr;

        const labelList& own = fvMesh_.faceOwner();
        const vectorField& Sf = fvMesh_.faceAreas();

        forAll(Sf, facei)
        {
            if (!fvMesh_.isInternalFace(facei))
            {
                magSb[own[facei]] += mag(Sf[facei]);
            }
        }
    }

    /*
    {
        scalarField& magSb = *magSbPtr;
        const polyBoundaryMesh& patches = fvMesh_.boundaryMesh();

        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];

            if (!(Pstream::parRun() && isA<processorPolyPatch>(pp)))
            {
                const labelUList& faceCells = pp.faceCells();
                const vectorField& pSf = pp.faceAreas();
                forAll(faceCells, pfi)
                {
                    magSb[faceCells[pfi]] += mag(pSf[pfi]);
                }
            }
        }
    }
    */

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
            *magSfPtr,
            *magSbPtr
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
            scalarField* aggMagSfPtr
            (
                new scalarField
                (
                    meshLevels_[nCreatedLevels].upperAddr().size(),
                    0
                )
            );

            restrictFaceField(*aggMagSfPtr, *magSfPtr, nCreatedLevels);

            if (nCreatedLevels)
            {
                delete magSfPtr;
            }

            magSfPtr = aggMagSfPtr;
        }

        // Agglomerate the cell boundary areas field for the next level
        {
            scalarField* aggMagSbPtr
            (
                new scalarField(meshLevels_[nCreatedLevels].size())
            );

            // Restrict but no parallel agglomeration (not supported)
            restrictField(*aggMagSbPtr, *magSbPtr, nCreatedLevels, false);

            delete magSbPtr;
            magSbPtr = aggMagSbPtr;
        }

        nCreatedLevels++;
    }

    // Shrink the storage of the levels to those created
    compactLevels(nCreatedLevels);

    // Delete temporary geometry storage
    if (nCreatedLevels)
    {
        delete VPtr;
        delete magSfPtr;
    }
    delete magSbPtr;
}


// ************************************************************************* //
