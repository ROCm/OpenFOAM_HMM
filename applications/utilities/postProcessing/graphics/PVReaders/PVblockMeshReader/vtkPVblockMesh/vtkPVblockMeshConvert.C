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

    const Map<string> blockStatus = getSelectedArrayMap
    (
        reader_->GetBlockSelection()
    );

    arrayRange& range = rangeBlocks_;
    range.block(blockNo);   // set output block
    label datasetNo = 0;    // restart at dataset 0

    const blockMesh& blkMesh = *meshPtr_;
    const pointField blkPoints(blkMesh.vertices() * blkMesh.scaleFactor());

    vtkIdType nodeIds[8];  // Space for VTK_HEXAHEDRON vertices
    int blockId = -1;
    for (auto partId : range)
    {
        ++blockId; // Increment first
        if (!blockStatus.found(partId))
        {
            continue;
        }
        const auto& longName = blockStatus[partId];

        const blockDescriptor& blockDef = blkMesh[blockId];
        const labelList& blockLabels = blockDef.blockShape();

        auto vtkpoints = vtkSmartPointer<vtkPoints>::New();
        vtkpoints->SetNumberOfPoints(blockLabels.size());

        forAll(blockLabels, pointi)
        {
            vtkpoints->SetPoint
            (
                pointi,
                blkPoints[blockLabels[pointi]].v_
            );
            nodeIds[pointi] = pointi;
        }

        auto vtkmesh = vtkSmartPointer<vtkUnstructuredGrid>::New();

        vtkmesh->Allocate(1);
        vtkmesh->InsertNextCell
        (
            VTK_HEXAHEDRON,
            8,
            nodeIds
        );

        vtkmesh->SetPoints(vtkpoints);

        addToBlock(output, vtkmesh, range, datasetNo, longName);
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
    const Map<string> edgeStatus = getSelectedArrayMap
    (
        reader_->GetCurvedEdgesSelection()
    );

    arrayRange& range = rangeEdges_;

    range.block(blockNo);      // set output block
    label datasetNo = 0;       // restart at dataset 0

    const blockMesh& blkMesh = *meshPtr_;
    const blockEdgeList& edges = blkMesh.edges();
    const scalar scaleFactor = blkMesh.scaleFactor();

    int edgeId = -1;
    for (auto partId : range)
    {
        ++edgeId; // Increment first
        if (!edgeStatus.found(partId))
        {
            continue;
        }
        const auto& longName = edgeStatus[partId];

        // Search each block
        forAll(blkMesh, blockId)
        {
            const blockDescriptor& blockDef = blkMesh[blockId];

            edgeList blkEdges = blockDef.blockShape().edges();

            // List of edge point and weighting factors
            pointField edgesPoints[12];
            scalarList edgesWeights[12];
            blockDef.edgesPointsWeights(edgesPoints, edgesWeights);

            // find the corresponding edge within the block
            label foundEdgeI = -1;
            forAll(blkEdges, blkEdgeI)
            {
                if (edges[edgeId].compare(blkEdges[blkEdgeI]))
                {
                    foundEdgeI = blkEdgeI;
                    break;
                }
            }

            if (foundEdgeI != -1)
            {
                const List<point>& edgePoints = edgesPoints[foundEdgeI];

                auto vtkpoints = vtkSmartPointer<vtkPoints>::New();
                vtkpoints->SetNumberOfPoints(edgePoints.size());

                vtkIdType pointIds[edgePoints.size()];
                forAll(edgePoints, pointi)
                {
                    const point p = edgePoints[pointi] * scaleFactor;

                    vtkpoints->SetPoint(pointi, p.v_);
                    pointIds[pointi] = pointi;
                }

                auto vtkmesh = vtkSmartPointer<vtkPolyData>::New();

                vtkmesh->Allocate(1);
                vtkmesh->InsertNextCell
                (
                    VTK_POLY_LINE,
                    edgePoints.size(),
                    pointIds
                );

                vtkmesh->SetPoints(vtkpoints);

                addToBlock(output, vtkmesh, range, datasetNo, longName);
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
        Info<< "<beg> " << FUNCTION_NAME << endl;
    }

    if (true)  // Or some flag or other condition
    {
        auto vtkpoints = vtkSmartPointer<vtkPoints>::New();
        vtkpoints->SetNumberOfPoints(blkPoints.size());

        forAll(blkPoints, pointi)
        {
            vtkpoints->SetPoint(pointi, blkPoints[pointi].v_);
        }

        auto vtkmesh = vtkSmartPointer<vtkPolyData>::New();

        vtkmesh->SetPoints(vtkpoints);
        vtkmesh->SetVerts(foamPvCore::identityVertices(blkPoints.size()));

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
        Info<< "<end> " << FUNCTION_NAME << endl;
    }
}


// ************************************************************************* //
