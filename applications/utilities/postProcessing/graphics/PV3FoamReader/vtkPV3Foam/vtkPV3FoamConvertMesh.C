/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
#include "fvMeshSubset.H"
#include "pointSet.H"
#include "vtkPV3FoamReader.h"

// VTK includes
#include "vtkDataArraySelection.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkUnstructuredGrid.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPV3Foam::convertMeshVolume
(
    vtkMultiBlockDataSet* output
)
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::convertMeshVolume" << endl;
    }

    // Create the internal mesh and add as DataSet(VOLUME, 0)
    if (selectedRegions_[VOLUME])
    {
        if (debug)
        {
            Info<< "Creating VTK internal mesh" << endl;
        }

        const fvMesh& mesh = *meshPtr_;

        vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::New();
        addVolumeMesh(mesh, ugrid, superCells_);
        AddToBlock(output, VOLUME, 0, ugrid, "internalMesh");
        SetBlockName(output, VOLUME, "Volume");
        selectedRegionDatasetIds_[VOLUME] = 0;
        ugrid->Delete();
    }
}


void Foam::vtkPV3Foam::convertMeshLagrangian
(
    vtkMultiBlockDataSet* output
)
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::convertMeshLagrangian" << endl;
    }

    // Create the Lagrangian mesh and add as DataSet(LAGRANGIAN, 0)
    if (lagrangianDataSize_)
    {
        if (selectedRegions_[LAGRANGIAN])
        {
            if (debug)
            {
                Info<< "Creating VTK Lagrangian mesh" << endl;
            }

            const fvMesh& mesh = *meshPtr_;

            vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::New();
            addLagrangianMesh(mesh, ugrid);
            AddToBlock(output, LAGRANGIAN, 0, ugrid, cloudName_.c_str());
            SetBlockName(output, LAGRANGIAN, "Lagrangian");
            selectedRegionDatasetIds_[LAGRANGIAN] = 0;
            ugrid->Delete();
        }
    }
}


void Foam::vtkPV3Foam::convertMeshPatches
(
    vtkMultiBlockDataSet* output
)
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::convertMeshPatches" << endl;
    }

    // Convert patches
    if (patchDataSize_)
    {
        const fvMesh& mesh = *meshPtr_;
        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        const label regionStartId = idRegionPatches_;
        const label regionEndId = idRegionPatches_ + patchDataSize_ - 1;

        // Create the patches and add as DataSet(VOLUME, ...)
        for (int i=regionStartId; i<=regionEndId; i++)
        {
            if (selectedRegions_[i])
            {
                const word regionName = reader_->GetRegionSelection()
                    ->GetArrayName(i);

                if (debug)
                {
                    Info<< "Creating VTK mesh for patch: " << regionName
                        << " region index: " << i << endl;
                }

                vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::New();
                const label patchId = mesh.boundaryMesh()
                    .findPatchID(regionName);
                addPatchMesh(patches[patchId], ugrid);
                const label nextId = GetNumberOfDataSets(output, VOLUME);
                AddToBlock(output, VOLUME, nextId, ugrid, regionName.c_str());
                selectedRegionDatasetIds_[i] = nextId;
                ugrid->Delete();
            }
        }
    }
}


void Foam::vtkPV3Foam::convertMeshCellSet
(
    vtkMultiBlockDataSet* output
)
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::convertMeshCellSet" << endl;
    }

    // Create the cell sets and add as DataSet(CELLSET, 0..n)
    if (cellSetDataSize_)
    {
        const fvMesh& mesh = *meshPtr_;

        const label regionStartId = idRegionCellSets_;
        const label regionEndId = regionStartId + cellSetDataSize_ - 1;

        for (int i=regionStartId; i<=regionEndId; i++)
        {
            if (selectedRegions_[i])
            {
                const word cSetName = reader_->GetRegionSelection()
                    ->GetArrayName(i);

                if (debug)
                {
                    Info<< "Creating VTK mesh for cellSet: " << cSetName
                        << " region index: " << i << endl;
                }

                const cellSet cSet(mesh, cSetName);
                fvMeshSubset subsetter(mesh);
                subsetter.setLargeCellSubset(cSet);

                vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::New();
                const label nextId = GetNumberOfDataSets(output, CELLSET);
                addVolumeMesh
                (
                    subsetter.subMesh(),
                    ugrid,
                    superCellSetCells_[nextId]
                );
                AddToBlock(output, CELLSET, nextId, ugrid, cSetName.c_str());
                selectedRegionDatasetIds_[i] = nextId;
                ugrid->Delete();
            }
        }

        SetBlockName(output, CELLSET, "CellSets");
    }
}


void Foam::vtkPV3Foam::convertMeshFaceSet
(
    vtkMultiBlockDataSet* output
)
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::convertMeshFaceSet" << endl;
    }

    // Create the face sets and add as DataSet(FACESET, 0..n)
    if (faceSetDataSize_)
    {
        const fvMesh& mesh = *meshPtr_;

        const label regionStartId = idRegionFaceSets_;
        const label regionEndId = regionStartId + faceSetDataSize_ - 1;

        for (int i=regionStartId; i<=regionEndId ; i++)
        {
            if (selectedRegions_[i])
            {
                const word fSetName = reader_->GetRegionSelection()
                    ->GetArrayName(i);

                if (debug)
                {
                    Info<< "Creating VTK mesh for faceSet: " << fSetName
                        << " region index: " << i << endl;
                }

                const faceSet fSet(mesh, fSetName);

                vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::New();
                addFaceSetMesh
                (
                    mesh,
                    fSet,
                    ugrid
                );
                const label nextId = GetNumberOfDataSets(output, FACESET);
                AddToBlock(output, FACESET, nextId, ugrid, fSetName.c_str());
                selectedRegionDatasetIds_[i] = nextId;
                ugrid->Delete();
            }
        }

        SetBlockName(output, FACESET, "FaceSets");
    }
}


void Foam::vtkPV3Foam::convertMeshPointSet
(
    vtkMultiBlockDataSet* output
)
{
    if (debug)
    {
        Info<< "entered Foam::vtkPV3Foam::convertMeshPointSet" << endl;
    }

    // Create the point sets and add as DataSet(POINTSET, 0..n)
    if (pointSetDataSize_)
    {
        const fvMesh& mesh = *meshPtr_;

        const label regionStartId = idRegionPointSets_;
        const label regionEndId = regionStartId + pointSetDataSize_ - 1;

        for (int i=regionStartId; i<=regionEndId ; i++)
        {
            if (selectedRegions_[i])
            {
                const word pSetName = reader_->GetRegionSelection()
                    ->GetArrayName(i);

                if (debug)
                {
                    Info<< "Creating VTK mesh for pointSet: " << pSetName
                        << " region index: " << i << endl;
                }

                const pointSet pSet(mesh, pSetName);

                vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::New();
                addPointSetMesh
                (
                    mesh,
                    pSet,
                    ugrid
                );
                label nextId = GetNumberOfDataSets(output, POINTSET);
                AddToBlock(output, POINTSET, nextId, ugrid, pSetName.c_str());
                selectedRegionDatasetIds_[i] = nextId;
                ugrid->Delete();
            }
        }

        SetBlockName(output, POINTSET, "PointSets");
    }
}


// ************************************************************************* //
