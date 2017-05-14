/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "vtkPVblockMesh.H"
#include "vtkPVblockMeshReader.h"

// OpenFOAM includes
#include "blockMesh.H"
#include "Time.H"

// VTK includes
#include "vtkCellArray.h"
#include "vtkDataArraySelection.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkSmartPointer.h"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPVblockMesh::convertMeshBlocks
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    if (debug)
    {
        Info<< "<beg> convertMeshBlocks" << endl;
    }

    vtkDataArraySelection* selection = reader_->GetBlockSelection();
    arrayRange& range = rangeBlocks_;
    range.block(blockNo);   // set output block
    label datasetNo = 0;    // restart at dataset 0

    const blockMesh& blkMesh = *meshPtr_;
    const pointField blkPoints(blkMesh.vertices() * blkMesh.scaleFactor());

    int blockI = 0;
    for
    (
        auto iter = range.cbegin();
        iter != range.cend();
        ++iter, ++blockI
    )
    {
        const auto partId = *iter;
        if (!blockStatus_[partId])
        {
            continue;
        }

        const blockDescriptor& blockDef = blkMesh[blockI];
        const labelList& blockLabels = blockDef.blockShape();

        vtkSmartPointer<vtkPoints> vtkpoints =
            vtkSmartPointer<vtkPoints>::New();

        vtkpoints->SetNumberOfPoints(blockLabels.size());

        vtkIdType nodeIds[8];
        forAll(blockLabels, pointi)
        {
            vtkpoints->SetPoint
            (
                pointi,
                blkPoints[blockLabels[pointi]].v_
            );
            nodeIds[pointi] = pointi;
        }

        vtkSmartPointer<vtkUnstructuredGrid> vtkmesh =
            vtkSmartPointer<vtkUnstructuredGrid>::New();

        vtkmesh->Allocate(1);
        vtkmesh->InsertNextCell
        (
            VTK_HEXAHEDRON,
            8,
            nodeIds
        );

        vtkmesh->SetPoints(vtkpoints);

        addToBlock
        (
            output, vtkmesh, range, datasetNo,
            selection->GetArrayName(partId)
        );
        ++datasetNo;
    }


    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> convertMeshBlocks" << endl;
    }
}


void Foam::vtkPVblockMesh::convertMeshEdges
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    vtkDataArraySelection* selection = reader_->GetCurvedEdgesSelection();
    arrayRange& range = rangeEdges_;

    range.block(blockNo);      // set output block
    label datasetNo = 0;       // restart at dataset 0

    const blockMesh& blkMesh = *meshPtr_;
    const blockEdgeList& edges = blkMesh.edges();
    const scalar scaleFactor = blkMesh.scaleFactor();

    int edgeI = 0;
    for
    (
        auto iter = range.cbegin();
        iter != range.cend();
        ++iter, ++edgeI
    )
    {
        const auto partId = *iter;
        if (!edgeStatus_[partId])
        {
            continue;
        }

        // search each block
        forAll(blkMesh, blockI)
        {
            const blockDescriptor& blockDef = blkMesh[blockI];

            edgeList blkEdges = blockDef.blockShape().edges();

            // List of edge point and weighting factors
            pointField edgesPoints[12];
            scalarList edgesWeights[12];
            blockDef.edgesPointsWeights(edgesPoints, edgesWeights);

            // find the corresponding edge within the block
            label foundEdgeI = -1;
            forAll(blkEdges, blkEdgeI)
            {
                if (edges[edgeI].compare(blkEdges[blkEdgeI]))
                {
                    foundEdgeI = blkEdgeI;
                    break;
                }
            }

            if (foundEdgeI != -1)
            {
                const List<point>& edgePoints = edgesPoints[foundEdgeI];

                vtkSmartPointer<vtkPoints> vtkpoints =
                    vtkSmartPointer<vtkPoints>::New();

                vtkpoints->SetNumberOfPoints(edgePoints.size());

                vtkIdType pointIds[edgePoints.size()];
                forAll(edgePoints, pointi)
                {
                    const point p = edgePoints[pointi] * scaleFactor;

                    vtkpoints->SetPoint(pointi, p.v_);
                    pointIds[pointi] = pointi;
                }

                vtkSmartPointer<vtkPolyData> vtkmesh =
                    vtkSmartPointer<vtkPolyData>::New();

                vtkmesh->Allocate(1);
                vtkmesh->InsertNextCell
                (
                    VTK_POLY_LINE,
                    edgePoints.size(),
                    pointIds
                );

                vtkmesh->SetPoints(vtkpoints);

                addToBlock
                (
                    output, vtkmesh, range, datasetNo,
                    selection->GetArrayName(partId)
                );
                ++datasetNo;

                break;
            }
        }
    }

    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> convertMeshEdges" << endl;
    }

}


void Foam::vtkPVblockMesh::convertMeshCorners
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = rangeCorners_;
    range.block(blockNo);      // set output block
    label datasetNo = 0;       // restart at dataset 0

    const pointField blkPoints(meshPtr_->vertices() * meshPtr_->scaleFactor());

    if (debug)
    {
        Info<< "<beg> convertMeshCorners" << endl;
    }

    if (true)  // or some flag or other condition
    {
        vtkSmartPointer<vtkPolyData> vtkmesh =
            vtkSmartPointer<vtkPolyData>::New();

        vtkSmartPointer<vtkPoints> vtkpoints =
            vtkSmartPointer<vtkPoints>::New();

        vtkSmartPointer<vtkCellArray> vtkcells =
            vtkSmartPointer<vtkCellArray>::New();

        vtkpoints->SetNumberOfPoints(blkPoints.size());
        vtkcells->Allocate(2*blkPoints.size());
        // If reusing memory, ensure insert always starts from 0
        vtkcells->Reset();

        vtkIdType pointId = 0;
        forAll(blkPoints, pointi)
        {
            vtkpoints->SetPoint(pointi, blkPoints[pointi].v_);
            vtkcells->InsertNextCell(1, &pointId);  // VTK_VERTEX
            pointId++;
        }

        vtkmesh->SetPoints(vtkpoints);
        vtkmesh->SetVerts(vtkcells);

        addToBlock(output, vtkmesh, range, datasetNo, range.name());
        ++datasetNo;
    }

    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> convertMeshCorners" << endl;
    }
}


// ************************************************************************* //
