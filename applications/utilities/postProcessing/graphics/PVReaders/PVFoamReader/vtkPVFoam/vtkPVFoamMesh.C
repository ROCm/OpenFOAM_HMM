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

#include "vtkPVFoam.H"

// OpenFOAM includes
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "fvMeshSubset.H"
#include "vtkPVFoamReader.h"
#include "uindirectPrimitivePatch.H"

// VTK includes
#include "vtkDataArraySelection.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPVFoam::convertMeshVolume
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = arrayRangeVolume_;
    range.block(blockNo);      // set output block
    label datasetNo = 0;       // restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    // resize for decomposed polyhedra
    regionPolyDecomp_.setSize(range.size());

    if (debug)
    {
        Info<< "<beg> convertMeshVolume" << endl;
        printMemory();
    }

    // Convert the internalMesh
    // this looks like more than one part, but it isn't
    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        const word partName = "internalMesh";

        if (!partStatus_[partId])
        {
            continue;
        }

        vtkUnstructuredGrid* vtkmesh = volumeVTKMesh
        (
            mesh,
            regionPolyDecomp_[datasetNo]
        );

        if (vtkmesh)
        {
            addToBlock(output, vtkmesh, range, datasetNo, partName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> convertMeshVolume" << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshLagrangian
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = arrayRangeLagrangian_;
    range.block(blockNo);      // set output block
    label datasetNo = 0;       // restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    if (debug)
    {
        Info<< "<beg> convertMeshLagrangian" << endl;
        printMemory();
    }

    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        const word cloudName = getPartName(partId);

        if (!partStatus_[partId])
        {
            continue;
        }

        vtkPolyData* vtkmesh = lagrangianVTKMesh(mesh, cloudName);

        if (vtkmesh)
        {
            addToBlock(output, vtkmesh, range, datasetNo, cloudName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> convertMeshLagrangian" << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshPatches
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = arrayRangePatches_;
    range.block(blockNo);      // set output block
    label datasetNo = 0;       // restart at dataset 0
    const fvMesh& mesh = *meshPtr_;
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (debug)
    {
        Info<< "<beg> convertMeshPatches" << endl;
        printMemory();
    }

    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        if (!partStatus_[partId])
        {
            continue;
        }

        const word patchName = getPartName(partId);

        labelHashSet
            patchIds(patches.patchSet(List<wordRe>(1, wordRe(patchName))));

        if (debug)
        {
            Info<< "Creating VTK mesh for patches [" << patchIds <<"] "
                << patchName << endl;
        }

        vtkPolyData* vtkmesh = nullptr;
        if (patchIds.size() == 1)
        {
            vtkmesh = patchVTKMesh(patchName, patches[patchIds.begin().key()]);
        }
        else
        {
            // Patch group. Collect patch faces.
            label sz = 0;
            forAllConstIter(labelHashSet, patchIds, iter)
            {
                sz += patches[iter.key()].size();
            }
            labelList faceLabels(sz);
            sz = 0;
            forAllConstIter(labelHashSet, patchIds, iter)
            {
                const polyPatch& pp = patches[iter.key()];
                forAll(pp, i)
                {
                    faceLabels[sz++] = pp.start()+i;
                }
            }

            uindirectPrimitivePatch pp
            (
                UIndirectList<face>
                (
                    mesh.faces(),
                    faceLabels
                ),
                mesh.points()
            );

            vtkmesh = patchVTKMesh(patchName, pp);
        }

        if (vtkmesh)
        {
            addToBlock(output, vtkmesh, range, datasetNo, patchName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> convertMeshPatches" << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshCellZones
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = arrayRangeCellZones_;
    range.block(blockNo);      // set output block
    label datasetNo = 0;       // restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    // resize for decomposed polyhedra
    zonePolyDecomp_.setSize(range.size());

    if (range.empty())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> convertMeshCellZones" << endl;
        printMemory();
    }

    const cellZoneMesh& zMesh = mesh.cellZones();
    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        const word zoneName = getPartName(partId);
        const label  zoneId = zMesh.findZoneID(zoneName);

        if (!partStatus_[partId] || zoneId < 0)
        {
            continue;
        }

        if (debug)
        {
            Info<< "Creating VTK mesh for cellZone[" << zoneId << "] "
                << zoneName << endl;
        }

        fvMeshSubset subsetter(mesh);
        subsetter.setLargeCellSubset(zMesh[zoneId]);

        vtkUnstructuredGrid* vtkmesh = volumeVTKMesh
        (
            subsetter.subMesh(),
            zonePolyDecomp_[datasetNo]
        );

        if (vtkmesh)
        {
            // superCells + addPointCellLabels must contain global cell ids
            inplaceRenumber
            (
                subsetter.cellMap(),
                zonePolyDecomp_[datasetNo].superCells()
            );
            inplaceRenumber
            (
                subsetter.cellMap(),
                zonePolyDecomp_[datasetNo].addPointCellLabels()
            );

            // copy pointMap as well, otherwise pointFields fail
            zonePolyDecomp_[datasetNo].pointMap() = subsetter.pointMap();

            addToBlock(output, vtkmesh, range, datasetNo, zoneName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> convertMeshCellZones" << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshCellSets
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = arrayRangeCellSets_;
    range.block(blockNo);      // set output block
    label datasetNo = 0;       // restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    // resize for decomposed polyhedra
    csetPolyDecomp_.setSize(range.size());

    if (debug)
    {
        Info<< "<beg> convertMeshCellSets" << endl;
        printMemory();
    }

    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        const word partName = getPartName(partId);

        if (!partStatus_[partId])
        {
            continue;
        }

        if (debug)
        {
            Info<< "Creating VTK mesh for cellSet=" << partName << endl;
        }

        const cellSet cSet(mesh, partName);
        fvMeshSubset subsetter(mesh);
        subsetter.setLargeCellSubset(cSet);

        vtkUnstructuredGrid* vtkmesh = volumeVTKMesh
        (
            subsetter.subMesh(),
            csetPolyDecomp_[datasetNo]
        );

        if (vtkmesh)
        {
            // superCells + addPointCellLabels must contain global cell ids
            inplaceRenumber
            (
                subsetter.cellMap(),
                csetPolyDecomp_[datasetNo].superCells()
            );
            inplaceRenumber
            (
                subsetter.cellMap(),
                csetPolyDecomp_[datasetNo].addPointCellLabels()
            );

            // copy pointMap as well, otherwise pointFields fail
            csetPolyDecomp_[datasetNo].pointMap() = subsetter.pointMap();

            addToBlock(output, vtkmesh, range, datasetNo, partName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> convertMeshCellSets" << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshFaceZones
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = arrayRangeFaceZones_;
    range.block(blockNo);      // set output block
    label datasetNo = 0;       // restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    if (range.empty())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> convertMeshFaceZones" << endl;
        printMemory();
    }

    const faceZoneMesh& zMesh = mesh.faceZones();
    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        const word zoneName = getPartName(partId);
        const label  zoneId = zMesh.findZoneID(zoneName);

        if (!partStatus_[partId] || zoneId < 0)
        {
            continue;
        }

        if (debug)
        {
            Info<< "Creating VTKmesh for faceZone[" << zoneId << "] "
                << zoneName << endl;
        }

        vtkPolyData* vtkmesh = patchVTKMesh(zoneName, zMesh[zoneId]());

        if (vtkmesh)
        {
            addToBlock(output, vtkmesh, range, datasetNo, zoneName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> convertMeshFaceZones" << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshFaceSets
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = arrayRangeFaceSets_;
    range.block(blockNo);      // set output block
    label datasetNo = 0;       // restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    if (debug)
    {
        Info<< "<beg> convertMeshFaceSets" << endl;
        printMemory();
    }

    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        const word partName = getPartName(partId);

        if (!partStatus_[partId])
        {
            continue;
        }

        if (debug)
        {
            Info<< "Creating VTK mesh for faceSet=" << partName << endl;
        }

        // faces in sorted order for more reliability
        uindirectPrimitivePatch p
        (
            UIndirectList<face>
            (
                mesh.faces(),
                faceSet(mesh, partName).sortedToc()
            ),
            mesh.points()
        );

        if (p.empty())
        {
            continue;
        }

        vtkPolyData* vtkmesh = patchVTKMesh("faceSet:" + partName, p);
        if (vtkmesh)
        {
            addToBlock(output, vtkmesh, range, datasetNo, partName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> convertMeshFaceSets" << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshPointZones
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = arrayRangePointZones_;
    range.block(blockNo);      // set output block
    label datasetNo = 0;       // restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    if (debug)
    {
        Info<< "<beg> convertMeshPointZones" << endl;
        printMemory();
    }

    if (range.size())
    {
        const pointZoneMesh& zMesh = mesh.pointZones();
        for (int partId = range.start(); partId < range.end(); ++partId)
        {
            word zoneName = getPartName(partId);
            label zoneId = zMesh.findZoneID(zoneName);

            if (!partStatus_[partId] || zoneId < 0)
            {
                continue;
            }

            const labelUList& pointLabels = zMesh[zoneId];

            vtkPoints* vtkpoints = vtkPoints::New();
            vtkpoints->Allocate(pointLabels.size());

            const pointField& meshPoints = mesh.points();
            forAll(pointLabels, pointi)
            {
                vtkpoints->InsertNextPoint(meshPoints[pointLabels[pointi]].v_);
            }

            vtkPolyData* vtkmesh = vtkPolyData::New();
            vtkmesh->SetPoints(vtkpoints);
            vtkpoints->Delete();

            if (vtkmesh)
            {
                addToBlock(output, vtkmesh, range, datasetNo, zoneName);
                vtkmesh->Delete();

                partDataset_[partId] = datasetNo++;
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
        Info<< "<end> convertMeshPointZones" << endl;
        printMemory();
    }
}



void Foam::vtkPVFoam::convertMeshPointSets
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = arrayRangePointSets_;
    range.block(blockNo);      // set output block
    label datasetNo = 0;       // restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    if (debug)
    {
        Info<< "<beg> convertMeshPointSets" << endl;
        printMemory();
    }

    for (int partId = range.start(); partId < range.end(); ++partId)
    {
        word partName = getPartName(partId);

        if (!partStatus_[partId])
        {
            continue;
        }

        if (debug)
        {
            Info<< "Creating VTK mesh for pointSet=" << partName << endl;
        }

        const pointSet pSet(mesh, partName);

        vtkPoints* vtkpoints = vtkPoints::New();
        vtkpoints->Allocate(pSet.size());

        const pointField& meshPoints = mesh.points();
        forAllConstIter(pointSet, pSet, iter)
        {
            vtkpoints->InsertNextPoint(meshPoints[iter.key()].v_);
        }

        vtkPolyData* vtkmesh = vtkPolyData::New();
        vtkmesh->SetPoints(vtkpoints);
        vtkpoints->Delete();

        if (vtkmesh)
        {
            addToBlock(output, vtkmesh, range, datasetNo, partName);
            vtkmesh->Delete();

            partDataset_[partId] = datasetNo++;
        }
    }

    // anything added?
    if (datasetNo)
    {
        ++blockNo;
    }

    if (debug)
    {
        Info<< "<end> convertMeshPointSets" << endl;
        printMemory();
    }
}


// ************************************************************************* //
