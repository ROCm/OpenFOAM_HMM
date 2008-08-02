/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description

\*---------------------------------------------------------------------------*/

#include "vtkPV3Foam.H"

// Foam includes
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "fvMeshSubset.H"
#include "vtkPV3FoamReader.h"

// VTK includes
#include "vtkDataArraySelection.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPV3Foam::convertMeshVolume
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::convertMeshVolume" << endl;
        printMemory();
    }

    const fvMesh& mesh = *meshPtr_;
    selectionInfo& selector = regionInfoVolume_;

    // set output block, restart at dataset 0
    selector.block(blockNo);
    label datasetNo = 0;

    // Create the internalMesh
    // TODO: multiple regions
    for
    (
        int regionId = selector.start();
        regionId < selector.end();
        ++regionId
    )
    {
        if (!regionStatus_[regionId])
        {
            continue;
        }

        if (debug)
        {
            Info<< "Creating VTK internalMesh" << endl;
        }

        vtkUnstructuredGrid* vtkmesh = volumeVTKMesh
        (
            mesh,
            superCells_
        );

        if (vtkmesh)
        {
            AddToBlock(output, selector, datasetNo, vtkmesh, "internalMesh");
            vtkmesh->Delete();

            regionDataset_[regionId] = datasetNo++;
        }
    }

    // was anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::convertMeshVolume" << endl;
        printMemory();
    }
}


void Foam::vtkPV3Foam::convertMeshLagrangian
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::convertMeshLagrangian" << endl;
        printMemory();
    }

    const fvMesh& mesh = *meshPtr_;
    selectionInfo& selector = regionInfoLagrangian_;
    vtkDataArraySelection* regionSelection = reader_->GetRegionSelection();

    // set output block, restart at dataset 0
    selector.block(blockNo);
    label datasetNo = 0;

    // Create Lagrangian meshes
    for
    (
        int regionId = selector.start();
        regionId < selector.end();
        ++regionId
    )
    {
        if (!regionStatus_[regionId])
        {
            continue;
        }

        word cloudName = getFirstWord
        (
            regionSelection->GetArrayName(regionId)
        );

        vtkPolyData* vtkmesh = lagrangianVTKMesh(mesh, cloudName);
        if (vtkmesh)
        {
            AddToBlock(output, selector, datasetNo, vtkmesh, cloudName);
            vtkmesh->Delete();

            regionDataset_[regionId] = datasetNo++;
        }
    }

    // was anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::convertMeshLagrangian" << endl;
        printMemory();
    }
}


void Foam::vtkPV3Foam::convertMeshPatches
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::convertMeshPatches" << endl;
        printMemory();
    }

    const fvMesh& mesh = *meshPtr_;
    selectionInfo& selector = regionInfoPatches_;
    vtkDataArraySelection* regionSelection = reader_->GetRegionSelection();

    // set output block, restart at dataset 0
    selector.block(blockNo);
    label datasetNo = 0;

    if (selector.size())
    {
        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        for
        (
            int regionId = selector.start();
            regionId < selector.end();
            ++regionId
        )
        {
            word patchName = getFirstWord
            (
                regionSelection->GetArrayName(regionId)
            );

            label patchId = patches.findPatchID(patchName);

            if (!regionStatus_[regionId] || patchId < 0)
            {
                continue;
            }

            if (debug)
            {
                Info<< "Creating VTK mesh for patch: " << patchName
                    << " patch index: " << patchId << endl;
            }

            vtkPolyData* vtkmesh = patchVTKMesh(patches[patchId]);
            if (vtkmesh)
            {
                AddToBlock(output, selector, datasetNo, vtkmesh, patchName);
                vtkmesh->Delete();

                regionDataset_[regionId] = datasetNo++;
            }
        }
    }

    // was anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::convertMeshPatches" << endl;
        printMemory();
    }
}


void Foam::vtkPV3Foam::convertMeshCellZones
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::convertMeshCellZones" << endl;
        printMemory();
    }

    const fvMesh& mesh = *meshPtr_;
    selectionInfo& selector = regionInfoCellZones_;
    vtkDataArraySelection* regionSelection = reader_->GetRegionSelection();

    // set output block, restart at dataset 0
    selector.block(blockNo);
    label datasetNo = 0;

    if (selector.size())
    {
        const cellZoneMesh& zMesh = mesh.cellZones();

        for
        (
            int regionId = selector.start();
            regionId < selector.end();
            ++regionId
        )
        {
            word zoneName = getFirstWord
            (
                regionSelection->GetArrayName(regionId)
            );

            label zoneId = zMesh.findZoneID(zoneName);

            if (!regionStatus_[regionId] || zoneId < 0)
            {
                continue;
            }

            if (debug)
            {
                Info<< "Creating VTK mesh for cellZone: "
                    << zoneId << endl;
            }

            fvMeshSubset subsetter(mesh);
            subsetter.setLargeCellSubset(zMesh[zoneId]);

            vtkUnstructuredGrid* vtkmesh = volumeVTKMesh
            (
                subsetter.subMesh(),
                zoneSuperCells_[datasetNo]
            );

            if (vtkmesh)
            {
                // renumber - superCells must contain global cell ids
                inplaceRenumber
                (
                    subsetter.cellMap(),
                    zoneSuperCells_[datasetNo]
                );

                AddToBlock
                (
                    output, selector, datasetNo, vtkmesh,
                    zMesh[zoneId].name() + ":cellZone"
                );
                vtkmesh->Delete();

                regionDataset_[regionId] = datasetNo++;
            }
        }
    }

    // was anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::convertMeshCellZones" << endl;
        printMemory();
    }
}


void Foam::vtkPV3Foam::convertMeshCellSets
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::convertMeshCellSets" << endl;
        printMemory();
    }

    const fvMesh& mesh = *meshPtr_;
    selectionInfo& selector = regionInfoCellSets_;
    vtkDataArraySelection* regionSelection = reader_->GetRegionSelection();

    // set output block, restart at dataset 0
    selector.block(blockNo);
    label datasetNo = 0;

    // Create the cell sets and add as dataset
    if (selector.size())
    {
        for
        (
            int regionId = selector.start();
            regionId < selector.end();
            ++regionId)
        {
            if (!regionStatus_[regionId])
            {
                continue;
            }

            word selectedName = getFirstWord
            (
                regionSelection->GetArrayName(regionId)
            );

            if (debug)
            {
                Info<< "Creating VTK mesh for cellSet: " << selectedName
                    << " region index: " << regionId << endl;
            }

            const cellSet cSet(mesh, selectedName);
            fvMeshSubset subsetter(mesh);
            subsetter.setLargeCellSubset(cSet);

            vtkUnstructuredGrid* vtkmesh = volumeVTKMesh
            (
                subsetter.subMesh(),
                csetSuperCells_[datasetNo]
            );

            if (vtkmesh)
            {
                // renumber - superCells must contain global cell ids
                inplaceRenumber
                (
                    subsetter.cellMap(),
                    csetSuperCells_[datasetNo]
                );

                AddToBlock
                (
                    output, selector, datasetNo, vtkmesh,
                    selectedName + ":cellSet"
                );
                vtkmesh->Delete();

                regionDataset_[regionId] = datasetNo++;
            }
        }
    }

    // was anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::convertMeshCellSets" << endl;
        printMemory();
    }
}


void Foam::vtkPV3Foam::convertMeshFaceZones
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::convertMeshFaceZones" << endl;
        printMemory();
    }

    const fvMesh& mesh = *meshPtr_;
    selectionInfo& selector = regionInfoFaceZones_;
    vtkDataArraySelection* regionSelection = reader_->GetRegionSelection();

    // set output block, restart at dataset 0
    selector.block(blockNo);
    label datasetNo = 0;

    // Create the cell zone(s)
    if (selector.size())
    {
        const faceZoneMesh& zMesh = mesh.faceZones();

        for
        (
            int regionId = selector.start();
            regionId < selector.end();
            ++regionId
        )
        {
            word zoneName = getFirstWord
            (
                regionSelection->GetArrayName(regionId)
            );

            const label zoneId = zMesh.findZoneID(zoneName);

            if (!regionStatus_[regionId] || zoneId < 0)
            {
                continue;
            }

            if (debug)
            {
                Info<< "Creating VTK mesh for faceZone[" << zoneId
                    << "] " << zoneName << endl;
            }

            vtkPolyData* vtkmesh = faceZoneVTKMesh
            (
                mesh,
                zMesh[zoneId]
            );

            if (vtkmesh)
            {
                AddToBlock
                (
                    output, selector, datasetNo, vtkmesh,
                    zMesh[zoneId].name() + ":faceZone"
                );
                vtkmesh->Delete();

                regionDataset_[regionId] = datasetNo++;
            }
        }
    }

    // was anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::convertMeshFaceZones" << endl;
        printMemory();
    }
}


void Foam::vtkPV3Foam::convertMeshFaceSets
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::convertMeshFaceSets" << endl;
        printMemory();
    }

    const fvMesh& mesh = *meshPtr_;
    selectionInfo& selector = regionInfoFaceSets_;
    vtkDataArraySelection* regionSelection = reader_->GetRegionSelection();

    // set output block, restart at dataset 0
    selector.block(blockNo);
    label datasetNo = 0;

    // Create the face sets and add as dataset
    if (selector.size())
    {
        for
        (
            int regionId = selector.start();
            regionId < selector.end();
            ++regionId
        )
        {
            if (!regionStatus_[regionId])
            {
                continue;
            }

            word selectedName = getFirstWord
            (
                regionSelection->GetArrayName(regionId)
            );

            if (debug)
            {
                Info<< "Creating VTK mesh for faceSet: " << selectedName
                    << " region index: " << regionId << endl;
            }

            const faceSet fSet(mesh, selectedName);

            vtkPolyData* vtkmesh = faceSetVTKMesh
            (
                mesh,
                fSet
            );

            if (vtkmesh)
            {
                AddToBlock
                (
                    output, selector, datasetNo, vtkmesh,
                    selectedName + ":faceSet"
                );
                vtkmesh->Delete();

                regionDataset_[regionId] = datasetNo++;
            }
        }
    }

    // was anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::convertMeshFaceSets" << endl;
        printMemory();
    }
}


void Foam::vtkPV3Foam::convertMeshPointZones
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::convertMeshPointZones" << endl;
        printMemory();
    }

    const fvMesh& mesh = *meshPtr_;
    selectionInfo& selector = regionInfoPointZones_;
    vtkDataArraySelection* regionSelection = reader_->GetRegionSelection();

    // set output block, restart at dataset 0
    selector.block(blockNo);
    label datasetNo = 0;

    if (selector.size())
    {
        const pointZoneMesh& zMesh = mesh.pointZones();

        for
        (
            int regionId = selector.start();
            regionId < selector.end();
            ++regionId
        )
        {
            word zoneName = getFirstWord
            (
                regionSelection->GetArrayName(regionId)
            );

            label zoneId = zMesh.findZoneID(zoneName);

            if (!regionStatus_[regionId] || zoneId < 0)
            {
                continue;
            }

            vtkPolyData* vtkmesh = pointZoneVTKMesh(mesh, zMesh[zoneId]);

            if (vtkmesh)
            {
                AddToBlock
                (
                    output, selector, datasetNo, vtkmesh,
                    zMesh[zoneId].name() + ":pointZone"
                );
                vtkmesh->Delete();

                regionDataset_[regionId] = datasetNo++;
            }
        }
    }

    // was anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::convertMeshPointZones" << endl;
        printMemory();
    }
}



void Foam::vtkPV3Foam::convertMeshPointSets
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::convertMeshPointSets" << endl;
        printMemory();
    }

    const fvMesh& mesh = *meshPtr_;
    selectionInfo& selector = regionInfoPointSets_;
    vtkDataArraySelection* regionSelection = reader_->GetRegionSelection();

    // set output block, restart at dataset 0
    selector.block(blockNo);
    label datasetNo = 0;

    if (selector.size())
    {
        for
        (
            int regionId = selector.start();
            regionId < selector.end();
            ++regionId
        )
        {
            if (!regionStatus_[regionId])
            {
                continue;
            }

            word selectedName = getFirstWord
            (
                regionSelection->GetArrayName(regionId)
            );

            if (debug)
            {
                Info<< "Creating VTK mesh for pointSet: " << selectedName
                    << " region index: " << regionId << endl;
            }

            const pointSet pSet(mesh, selectedName);

            vtkPolyData* vtkmesh = pointSetVTKMesh(mesh, pSet);
            if (vtkmesh)
            {
                AddToBlock
                (
                    output, selector, datasetNo, vtkmesh,
                    selectedName + ":pointSet"
                );
                vtkmesh->Delete();

                regionDataset_[regionId] = datasetNo++;
            }
        }
    }

    // was anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::convertMeshPointSets" << endl;
        printMemory();
    }
}

// ************************************************************************* //
