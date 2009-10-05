/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "vtkPV3FoamBlockMesh.H"
#include "vtkPV3FoamBlockMeshReader.h"

// Foam includes
#include "blockMesh.H"
#include "Time.H"

#include "vtkOpenFOAMPoints.H"

// VTK includes
#include "vtkCellArray.h"
#include "vtkDataArraySelection.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPV3FoamBlockMesh::convertMeshBlocks
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    partInfo& selector = partInfoBlocks_;
    selector.block(blockNo);   // set output block
    label datasetNo = 0;       // restart at dataset 0
    const blockMesh& blkMesh = *meshPtr_;

    const Foam::pointField& blockPoints = blkMesh.blockPointField();


    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3FoamBlockMesh::convertMeshBlocks" << endl;
    }

    int blockI = 0;

    for
    (
        int partId = selector.start();
        partId < selector.end();
        ++partId, ++blockI
    )
    {
        if (!partStatus_[partId])
        {
            continue;
        }

        const word partName = Foam::name(blockI);
        const blockDescriptor& blockDef = blkMesh[blockI].blockDef();


        vtkUnstructuredGrid* vtkmesh = vtkUnstructuredGrid::New();

        // Convert Foam mesh vertices to VTK
        vtkPoints *vtkpoints = vtkPoints::New();
        vtkpoints->Allocate( blockDef.nPoints() );
        const labelList& blockLabels = blockDef.blockShape();

        vtkmesh->Allocate(1);
        vtkIdType nodeIds[8];

        forAll(blockLabels, ptI)
        {
            vtkInsertNextOpenFOAMPoint
            (
                vtkpoints,
                blockPoints[blockLabels[ptI]]
            );

            nodeIds[ptI] = ptI;
        }

        vtkmesh->InsertNextCell
        (
            VTK_HEXAHEDRON,
            8,
            nodeIds
        );

        vtkmesh->SetPoints(vtkpoints);
        vtkpoints->Delete();


        AddToBlock(output, vtkmesh, selector, datasetNo, partName);
        vtkmesh->Delete();

        partDataset_[partId] = datasetNo++;
    }


    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3FoamBlockMesh::convertMeshBlocks" << endl;
    }
}


void Foam::vtkPV3FoamBlockMesh::convertMeshCorners
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    partInfo& selector = partInfoCorners_;

    selector.block(blockNo);   // set output block
    label datasetNo = 0;       // restart at dataset 0

    const pointField& blockPoints = meshPtr_->blockPointField();

    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3FoamBlockMesh::convertMeshCorners" << endl;
    }

    if (true)  // or some flag or other condition
    {
        vtkPolyData* vtkmesh = vtkPolyData::New();
        vtkPoints* vtkpoints = vtkPoints::New();
        vtkCellArray* vtkcells = vtkCellArray::New();

        vtkpoints->Allocate( blockPoints.size() );
        vtkcells->Allocate( blockPoints.size() );

        vtkIdType pointId = 0;
        forAll(blockPoints, ptI)
        {
            vtkInsertNextOpenFOAMPoint(vtkpoints, blockPoints[ptI]);

            vtkcells->InsertNextCell(1, &pointId);
            pointId++;
        }

        vtkmesh->SetPoints(vtkpoints);
        vtkpoints->Delete();

        vtkmesh->SetVerts(vtkcells);
        vtkcells->Delete();

        AddToBlock(output, vtkmesh, selector, datasetNo, "");
        vtkmesh->Delete();

        datasetNo++;
    }

    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3FoamBlockMesh::convertMeshCorners" << endl;
    }
}


// ************************************************************************* //
