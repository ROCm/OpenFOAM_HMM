/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

#include "fvcSmooth.H"
#include "volFields.H"
#include "FaceCellWave.H"
#include "smoothData.H"
#include "sweepData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::scalar Foam::smoothData::maxRatio = 1.0;

void Foam::fvc::smooth
(
    volScalarField& field,
    const scalar coeff
)
{
    const fvMesh& mesh = field.mesh();
    smoothData::maxRatio = 1 + coeff;

    DynamicList<label> changedFaces(mesh.nFaces()/100 + 100);
    DynamicList<smoothData> changedFacesInfo(changedFaces.size());

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nbr = neighbour[facei];

        // Check if owner value much larger than neighbour value or vice versa
        if (field[own] > smoothData::maxRatio*field[nbr])
        {
            changedFaces.append(facei);
            changedFacesInfo.append(smoothData(field[own]));
        }
        else if (field[nbr] > smoothData::maxRatio*field[own])
        {
            changedFaces.append(facei);
            changedFacesInfo.append(smoothData(field[nbr]));
        }
    }

    // Insert all faces of coupled patches - FaceCellWave will correct them
    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        if (patch.coupled())
        {
            forAll(patch, patchFacei)
            {
                label facei = patch.start() + patchFacei;
                label own = mesh.faceOwner()[facei];

                changedFaces.append(facei);
                changedFacesInfo.append(smoothData(field[own]));
            }
        }
    }

    changedFaces.shrink();
    changedFacesInfo.shrink();

    // Set initial field on cells
    List<smoothData> cellData(mesh.nCells());

    forAll(field, celli)
    {
        cellData[celli] = field[celli];
    }

    // Set initial field on faces
    List<smoothData> faceData(mesh.nFaces());


    // Propagate information over whole domain
    FaceCellWave<smoothData > smoothData
    (
        mesh,
        changedFaces,
        changedFacesInfo,
        faceData,
        cellData,
        mesh.globalData().nTotalCells()    // max iterations
    );

    forAll(field, celli)
    {
        field[celli] = cellData[celli].value();
    }

    field.correctBoundaryConditions();
}


void Foam::fvc::spread
(
    volScalarField& field,
    const volScalarField& alpha,
    const label nLayers
)
{
    const fvMesh& mesh = field.mesh();
    smoothData::maxRatio = 1;

    DynamicList<label> changedFaces(mesh.nFaces()/100 + 100);
    DynamicList<smoothData> changedFacesInfo(changedFaces.size());

    // Set initial field on cells
    List<smoothData> cellData(mesh.nCells());

    forAll(field, celli)
    {
        cellData[celli] = field[celli];
    }

    // Set initial field on faces
    List<smoothData> faceData(mesh.nFaces());

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nbr = neighbour[facei];

        if
        (
            (alpha[own] > 0.01 && alpha[own] < 0.99)
         || (alpha[nbr] > 0.01 && alpha[nbr] < 0.99)
        )
        {
            if (mag(alpha[own] - alpha[nbr]) > 0.2)
            {
                changedFaces.append(facei);
                changedFacesInfo.append
                (
                    smoothData(max(field[own], field[nbr]))
                );
            }
        }
    }

    // Insert all faces of coupled patches - FaceCellWave will correct them
    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        if (patch.coupled())
        {
            forAll(patch, patchFacei)
            {
                label facei = patch.start() + patchFacei;
                label own = mesh.faceOwner()[facei];

                scalarField alphapn =
                    alpha.boundaryField()[patchi].patchNeighbourField();

                if
                (
                    (alpha[own] > 0.01 && alpha[own] < 0.99)
                 || (alphapn[patchFacei] > 0.01 && alphapn[patchFacei] < 0.99)
                )
                {
                    if (mag(alpha[own] - alphapn[patchFacei]) > 0.2)
                    {
                        changedFaces.append(facei);
                        changedFacesInfo.append(smoothData(field[own]));
                    }
                }
            }
        }
    }

    changedFaces.shrink();
    changedFacesInfo.shrink();

    // Propagate information over whole domain
    FaceCellWave<smoothData> smoothData
    (
        mesh,
        faceData,
        cellData
    );

    smoothData.setFaceInfo(changedFaces, changedFacesInfo);

    smoothData.iterate(nLayers);

    forAll(field, celli)
    {
        field[celli] = cellData[celli].value();
    }

    field.correctBoundaryConditions();
}


void Foam::fvc::sweep
(
    volScalarField& field,
    const volScalarField& alpha,
    const label nLayers
)
{
    const fvMesh& mesh = field.mesh();

    DynamicList<label> changedFaces(mesh.nFaces()/100 + 100);
    DynamicList<sweepData> changedFacesInfo(changedFaces.size());

    // Set initial field on cells
    List<sweepData> cellData(mesh.nCells());

    // Set initial field on faces
    List<sweepData> faceData(mesh.nFaces());

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();
    const vectorField& Cf = mesh.faceCentres();

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nbr = neighbour[facei];

        if (mag(alpha[own] - alpha[nbr]) > 0.2)
        {
            changedFaces.append(facei);
            changedFacesInfo.append
            (
                sweepData(max(field[own], field[nbr]), Cf[facei])
            );
        }
    }

    // Insert all faces of coupled patches - FaceCellWave will correct them
    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        if (patch.coupled())
        {
            forAll(patch, patchFacei)
            {
                label facei = patch.start() + patchFacei;
                label own = mesh.faceOwner()[facei];

                scalarField alphapn =
                    alpha.boundaryField()[patchi].patchNeighbourField();

                if (mag(alpha[own] - alphapn[patchFacei]) > 0.2)
                {
                    changedFaces.append(facei);
                    changedFacesInfo.append
                    (
                        sweepData(field[own], Cf[facei])
                    );
                }
            }
        }
    }

    changedFaces.shrink();
    changedFacesInfo.shrink();

    // Propagate information over whole domain
    FaceCellWave<sweepData> sweepData
    (
        mesh,
        faceData,
        cellData
    );

    sweepData.setFaceInfo(changedFaces, changedFacesInfo);

    sweepData.iterate(nLayers);

    forAll(field, celli)
    {
        if (cellData[celli].valid())
        {
            field[celli] = max(field[celli], cellData[celli].value());
        }
    }

    field.correctBoundaryConditions();
}


// ************************************************************************* //
