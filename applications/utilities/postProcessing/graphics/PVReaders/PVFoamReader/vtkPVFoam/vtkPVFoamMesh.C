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
#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPVFoam::convertMeshVolume
(
    vtkMultiBlockDataSet* output,
    int& blockNo
)
{
    arrayRange& range = rangeVolume_;
    range.block(blockNo);      // set output block
    label datasetNo = 0;       // restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    if (debug)
    {
        Info<< "<beg> convertMeshVolume" << endl;
        printMemory();
    }

    // Convert the internalMesh
    // this looks like more than one part, but it isn't
    for (auto partId : range)
    {
        if (!selectedPartIds_.found(partId))
        {
            continue;
        }

        const auto& longName = selectedPartIds_[partId];
        const word partName = getPartName(partId);

        vtkSmartPointer<vtkUnstructuredGrid> vtkmesh = volumeVTKMesh
        (
            mesh,
            cachedVtu_(longName)
        );

        if (vtkmesh)
        {
            addToBlock(output, vtkmesh, range, datasetNo, partName);
            partDataset_.set(partId, datasetNo++);
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
    arrayRange& range = rangeLagrangian_;
    range.block(blockNo);      // set output block
    label datasetNo = 0;       // restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    if (debug)
    {
        Info<< "<beg> convertMeshLagrangian" << endl;
        printMemory();
    }

    for (auto partId : range)
    {
        if (!selectedPartIds_.found(partId))
        {
            continue;
        }

        const word cloudName = getPartName(partId);

        vtkSmartPointer<vtkPolyData> vtkmesh =
            lagrangianVTKMesh(mesh, cloudName);

        if (vtkmesh)
        {
            addToBlock(output, vtkmesh, range, datasetNo, cloudName);
            partDataset_.set(partId, datasetNo++);
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
    arrayRange& range = rangePatches_;
    range.block(blockNo);      // set output block
    label datasetNo = 0;       // restart at dataset 0
    const fvMesh& mesh = *meshPtr_;
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (debug)
    {
        Info<< "<beg> convertMeshPatches" << endl;
        printMemory();
    }

    for (auto partId : range)
    {
        if (!selectedPartIds_.found(partId))
        {
            continue;
        }
        const auto& longName = selectedPartIds_[partId];
        const word partName = getPartName(partId);

        vtkSmartPointer<vtkPolyData> vtkmesh;

        if (longName.startsWith("group/"))
        {
            // Patch group. Collect patch faces.

            const labelList& patchIds =
                patches.groupPatchIDs().lookup(partName, labelList());

            if (debug)
            {
                Info<< "Creating VTK mesh for patches [" << patchIds <<"] "
                    << longName << endl;
            }

            label sz = 0;
            for (auto id : patchIds)
            {
                sz += patches[id].size();
            }
            labelList faceLabels(sz);
            sz = 0;
            for (auto id : patchIds)
            {
                const auto& pp = patches[id];
                forAll(pp, i)
                {
                    faceLabels[sz++] = pp.start()+i;
                }
            }

            if (faceLabels.size())
            {
                uindirectPrimitivePatch pp
                (
                    UIndirectList<face>(mesh.faces(), faceLabels),
                    mesh.points()
                );

                vtkmesh = patchVTKMesh(partName, pp);
            }
        }
        else
        {
            const label patchId = patches.findPatchID(partName);

            if (debug)
            {
                Info<< "Creating VTK mesh for patch [" << patchId <<"] "
                    << partName << endl;
            }

            if (patchId >= 0)
            {
                vtkmesh = patchVTKMesh(partName, patches[patchId]);
            }
        }

        if (vtkmesh)
        {
            addToBlock(output, vtkmesh, range, datasetNo, partName);
            partDataset_.set(partId, datasetNo++);
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
    arrayRange& range = rangeCellZones_;
    range.block(blockNo);      // set output block
    label datasetNo = 0;       // restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

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
    for (auto partId : range)
    {
        if (!selectedPartIds_.found(partId))
        {
            continue;
        }

        const auto& longName = selectedPartIds_[partId];
        const word zoneName = getPartName(partId);
        const label  zoneId = zMesh.findZoneID(zoneName);

        if (zoneId < 0)
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

        vtkSmartPointer<vtkUnstructuredGrid> vtkmesh = volumeVTKMesh
        (
            subsetter.subMesh(),
            cachedVtu_(longName)
        );

        if (vtkmesh)
        {
            // Convert cellMap, addPointCellLabels to global cell ids
            inplaceRenumber
            (
                subsetter.cellMap(),
                cachedVtu_[longName].cellMap()
            );
            inplaceRenumber
            (
                subsetter.cellMap(),
                cachedVtu_[longName].additionalIds()
            );

            // copy pointMap as well, otherwise pointFields fail
            cachedVtu_[longName].pointMap() = subsetter.pointMap();

            addToBlock(output, vtkmesh, range, datasetNo, zoneName);
            partDataset_.set(partId, datasetNo++);
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
    arrayRange& range = rangeCellSets_;
    range.block(blockNo);      // set output block
    label datasetNo = 0;       // restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    if (debug)
    {
        Info<< "<beg> convertMeshCellSets" << endl;
        printMemory();
    }

    for (auto partId : range)
    {
        if (!selectedPartIds_.found(partId))
        {
            continue;
        }

        const auto& longName = selectedPartIds_[partId];
        const word partName = getPartName(partId);

        if (debug)
        {
            Info<< "Creating VTK mesh for cellSet=" << partName << endl;
        }

        const cellSet cSet(mesh, partName);
        fvMeshSubset subsetter(mesh);
        subsetter.setLargeCellSubset(cSet);

        vtkSmartPointer<vtkUnstructuredGrid> vtkmesh = volumeVTKMesh
        (
            subsetter.subMesh(),
            cachedVtu_(longName)
        );

        if (vtkmesh)
        {
            // Convert cellMap, addPointCellLabels to global cell ids
            inplaceRenumber
            (
                subsetter.cellMap(),
                cachedVtu_[longName].cellMap()
            );
            inplaceRenumber
            (
                subsetter.cellMap(),
                cachedVtu_[longName].additionalIds()
            );

            // copy pointMap as well, otherwise pointFields fail
            cachedVtu_[longName].pointMap() = subsetter.pointMap();

            addToBlock(output, vtkmesh, range, datasetNo, partName);
            partDataset_.set(partId, datasetNo++);
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
    arrayRange& range = rangeFaceZones_;
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
    for (auto partId : range)
    {
        if (!selectedPartIds_.found(partId))
        {
            continue;
        }

        const word zoneName = getPartName(partId);
        const label  zoneId = zMesh.findZoneID(zoneName);

        if (zoneId < 0)
        {
            continue;
        }

        if (debug)
        {
            Info<< "Creating VTKmesh for faceZone[" << zoneId << "] "
                << zoneName << endl;
        }

        vtkSmartPointer<vtkPolyData> vtkmesh =
            patchVTKMesh(zoneName, zMesh[zoneId]());

        if (vtkmesh)
        {
            addToBlock(output, vtkmesh, range, datasetNo, zoneName);
            partDataset_.set(partId, datasetNo++);
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
    arrayRange& range = rangeFaceSets_;
    range.block(blockNo);      // set output block
    label datasetNo = 0;       // restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    if (debug)
    {
        Info<< "<beg> convertMeshFaceSets" << endl;
        printMemory();
    }

    for (auto partId : range)
    {
        if (!selectedPartIds_.found(partId))
        {
            continue;
        }

        const word partName = getPartName(partId);

        if (debug)
        {
            Info<< "Creating VTK mesh for faceSet=" << partName << endl;
        }

        // faces in sorted order for more reliability
        const labelList faceLabels = faceSet(mesh, partName).sortedToc();

        uindirectPrimitivePatch p
        (
            UIndirectList<face>(mesh.faces(), faceLabels),
            mesh.points()
        );

        if (p.empty())
        {
            continue;
        }

        vtkSmartPointer<vtkPolyData> vtkmesh =
            patchVTKMesh("faceSet:" + partName, p);

        if (vtkmesh)
        {
            addToBlock(output, vtkmesh, range, datasetNo, partName);
            partDataset_.set(partId, datasetNo++);
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
    arrayRange& range = rangePointZones_;
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
        for (auto partId : range)
        {
            if (!selectedPartIds_.found(partId))
            {
                continue;
            }

            const word zoneName = getPartName(partId);
            const label zoneId = zMesh.findZoneID(zoneName);

            if (zoneId < 0)
            {
                continue;
            }

            const pointField& meshPoints = mesh.points();
            const labelUList& pointLabels = zMesh[zoneId];

            vtkSmartPointer<vtkPoints> vtkpoints =
                vtkSmartPointer<vtkPoints>::New();

            vtkpoints->SetNumberOfPoints(pointLabels.size());

            forAll(pointLabels, pointi)
            {
                vtkpoints->SetPoint
                (
                    pointi,
                    meshPoints[pointLabels[pointi]].v_
                );
            }

            vtkSmartPointer<vtkPolyData> vtkmesh =
                vtkSmartPointer<vtkPolyData>::New();

            vtkmesh->SetPoints(vtkpoints);

            if (vtkmesh)
            {
                addToBlock(output, vtkmesh, range, datasetNo, zoneName);
                partDataset_.set(partId, datasetNo++);
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
    arrayRange& range = rangePointSets_;
    range.block(blockNo);      // set output block
    label datasetNo = 0;       // restart at dataset 0
    const fvMesh& mesh = *meshPtr_;

    if (debug)
    {
        Info<< "<beg> convertMeshPointSets" << endl;
        printMemory();
    }

    for (auto partId : range)
    {
        if (!selectedPartIds_.found(partId))
        {
            continue;
        }

        const word partName = getPartName(partId);

        if (debug)
        {
            Info<< "Creating VTK mesh for pointSet=" << partName << endl;
        }

        const pointField& meshPoints = mesh.points();
        const labelList pointLabels = pointSet(mesh, partName).sortedToc();

        vtkSmartPointer<vtkPoints> vtkpoints =
            vtkSmartPointer<vtkPoints>::New();

        vtkpoints->SetNumberOfPoints(pointLabels.size());

        forAll(pointLabels, pointi)
        {
            vtkpoints->SetPoint
            (
                pointi,
                meshPoints[pointLabels[pointi]].v_
            );
        }

        vtkSmartPointer<vtkPolyData> vtkmesh =
            vtkSmartPointer<vtkPolyData>::New();

        vtkmesh->SetPoints(vtkpoints);

        if (vtkmesh)
        {
            addToBlock(output, vtkmesh, range, datasetNo, partName);
            partDataset_.set(partId, datasetNo++);
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
