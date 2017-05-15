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

void Foam::vtkPVFoam::convertMeshVolume()
{
    arrayRange& range = rangeVolume_;
    const fvMesh& mesh = *meshPtr_;

    if (debug)
    {
        Info<< "<beg> " << FUNCTION_NAME << endl;
        printMemory();
    }

    // Convert the internalMesh
    // this looks like more than one part, but it isn't
    for (auto partId : range)
    {
        if (selectedPartIds_.found(partId))
        {
            const auto& longName = selectedPartIds_[partId];

            vtkSmartPointer<vtkUnstructuredGrid> vtkmesh =
                volumeVTKMesh
                (
                    mesh,
                    cachedVtu_(longName)
                );

            cachedVtu_[longName].vtkmesh = vtkmesh;
        }
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshLagrangian()
{
    arrayRange& range = rangeLagrangian_;
    const fvMesh& mesh = *meshPtr_;

    if (debug)
    {
        Info<< "<beg> " << FUNCTION_NAME << endl;
        printMemory();
    }

    for (auto partId : range)
    {
        if (selectedPartIds_.found(partId))
        {
            const auto& longName = selectedPartIds_[partId];
            const word cloudName = getPartName(partId);

            vtkSmartPointer<vtkPolyData> vtkmesh =
                lagrangianVTKMesh
                (
                    mesh,
                    cloudName
                );

            cachedVtp_(longName).vtkmesh = vtkmesh;
        }
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshPatches()
{
    arrayRange& range = rangePatches_;
    const fvMesh& mesh = *meshPtr_;
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (debug)
    {
        Info<< "<beg> " << FUNCTION_NAME << endl;
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
            cachedVtp_(longName).vtkmesh = vtkmesh;
        }
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshSubset
(
    const fvMeshSubset& subsetter,
    const string& longName
)
{
    vtkSmartPointer<vtkUnstructuredGrid> vtkmesh = volumeVTKMesh
    (
        subsetter.subMesh(),
        cachedVtu_(longName)
    );

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
    cachedVtu_[longName].vtkmesh = vtkmesh;
}


void Foam::vtkPVFoam::convertMeshCellZones()
{
    const arrayRange& range = rangeCellZones_;
    const fvMesh& mesh = *meshPtr_;

    if (range.empty())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> " << FUNCTION_NAME << endl;
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

        convertMeshSubset(subsetter, longName);
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshCellSets()
{
    const arrayRange& range = rangeCellSets_;
    const fvMesh& mesh = *meshPtr_;

    if (debug)
    {
        Info<< "<beg> " << FUNCTION_NAME << endl;
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

        convertMeshSubset(subsetter, longName);
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshFaceZones()
{
    arrayRange& range = rangeFaceZones_;
    const fvMesh& mesh = *meshPtr_;

    if (range.empty())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> " << FUNCTION_NAME << endl;
        printMemory();
    }

    const faceZoneMesh& zMesh = mesh.faceZones();
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
            Info<< "Creating VTKmesh for faceZone[" << zoneId << "] "
                << zoneName << endl;
        }

        vtkSmartPointer<vtkPolyData> vtkmesh =
            patchVTKMesh
            (
                zoneName,
                zMesh[zoneId]()
            );

        cachedVtp_(longName).vtkmesh = vtkmesh;
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshFaceSets()
{
    arrayRange& range = rangeFaceSets_;
    const fvMesh& mesh = *meshPtr_;

    if (debug)
    {
        Info<< "<beg> " << FUNCTION_NAME << endl;
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
            patchVTKMesh
            (
                "faceSet:" + partName,
                p
            );

        cachedVtp_(longName).vtkmesh = vtkmesh;
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshPointZones()
{
    arrayRange& range = rangePointZones_;
    const fvMesh& mesh = *meshPtr_;

    if (debug)
    {
        Info<< "<beg> " << FUNCTION_NAME << endl;
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

            const auto& longName = selectedPartIds_[partId];
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

            cachedVtp_(longName).vtkmesh = vtkmesh;
        }
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshPointSets()
{
    arrayRange& range = rangePointSets_;
    const fvMesh& mesh = *meshPtr_;

    if (debug)
    {
        Info<< "<beg> " << FUNCTION_NAME << endl;
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

        cachedVtp_(longName).vtkmesh = vtkmesh;
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << endl;
        printMemory();
    }
}


// ************************************************************************* //
