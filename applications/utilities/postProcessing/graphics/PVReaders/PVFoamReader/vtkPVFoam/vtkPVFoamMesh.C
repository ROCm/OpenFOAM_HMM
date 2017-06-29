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
    const arrayRange& range = rangeVolume_;
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
            foamVtuData& vtuData = cachedVtu_(longName);

            vtkSmartPointer<vtkUnstructuredGrid> vtkgeom;
            if (vtuData.nPoints())
            {
                if (meshState_ == polyMesh::UNCHANGED)
                {
                    if (debug)
                    {
                        Info<< "reuse " << longName << nl;
                    }
                    vtuData.reuse(); // reuse
                    continue;
                }
                else if (meshState_ == polyMesh::POINTS_MOVED)
                {
                    if (debug)
                    {
                        Info<< "move points " << longName << nl;
                    }
                    vtkgeom = vtuData.getCopy();
                    vtkgeom->SetPoints(movePoints(mesh, vtuData));
                }
            }

            if (!vtkgeom)
            {
                // Nothing usable from cache - create new geometry
                vtkgeom = volumeVTKMesh(mesh, vtuData);
            }

            vtuData.set(vtkgeom);
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
    const arrayRange& range = rangeLagrangian_;
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
            const word cloudName = getFoamName(longName);

            // Direct conversion, no caching for Lagrangian
            cachedVtp_(longName).set
            (
                lagrangianVTKMesh(mesh, cloudName)
            );
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
    const arrayRange& range = rangePatches_;
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
        const word partName = getFoamName(longName);
        foamVtpData& vtpData = cachedVtp_(longName);

        vtkSmartPointer<vtkPolyData> vtkgeom;
        if (vtpData.nPoints())
        {
            if (meshState_ == polyMesh::UNCHANGED)
            {
                // Without movement is easy.
                if (debug)
                {
                    Info<< "reuse " << longName << nl;
                }
                vtpData.reuse();
                continue;
            }
            else if (meshState_ == polyMesh::POINTS_MOVED)
            {
                // Point movement on single patch is OK

                const labelList& patchIds = vtpData.additionalIds();
                if (patchIds.size() == 1)
                {
                    vtkgeom = vtpData.getCopy();
                    vtkgeom->SetPoints(movePatchPoints(patches[patchIds[0]]));
                    continue;
                }
            }
        }

        if (longName.startsWith("group/"))
        {
            // Patch group. Collect patch faces.

            vtpData.clear(); // Remove any old mappings

            const labelList& patchIds =
                patches.groupPatchIDs().lookup(partName, labelList());

            if (debug)
            {
                Info<< "Creating VTK mesh for patches [" << patchIds <<"] "
                    << longName << endl;
            }

            // Store good patch ids as additionalIds
            vtpData.additionalIds().reserve(patchIds.size());

            label sz = 0;
            for (auto id : patchIds)
            {
                const label n = patches[id].size();

                if (n)
                {
                    vtpData.additionalIds().append(id);
                    sz += n;
                }
            }
            Foam::sort(vtpData.additionalIds());

            // Temporarily (mis)use cellMap for face labels
            DynamicList<label>& faceLabels = vtpData.cellMap();
            faceLabels.reserve(sz);

            for (auto id : vtpData.additionalIds())
            {
                const auto& pp = patches[id];
                forAll(pp, i)
                {
                    faceLabels.append(pp.start()+i);
                }
            }

            if (faceLabels.size())
            {
                uindirectPrimitivePatch pp
                (
                    UIndirectList<face>(mesh.faces(), faceLabels),
                    mesh.points()
                );

                vtkgeom = patchVTKMesh(pp);
            }

            faceLabels.clear();  // Unneeded
        }
        else
        {
            vtpData.clear(); // Remove any old mappings

            const label patchId = patches.findPatchID(partName);

            if (debug)
            {
                Info<< "Creating VTK mesh for patch [" << patchId <<"] "
                    << partName << endl;
            }

            if (patchId >= 0)
            {
                // Store good patch id as additionalIds
                vtpData.additionalIds() = {patchId};

                vtkgeom = patchVTKMesh(patches[patchId]);
            }
        }

        if (vtkgeom)
        {
            vtpData.set(vtkgeom);
        }
        else
        {
            cachedVtp_.erase(longName);
        }
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << endl;
        printMemory();
    }
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
        const word zoneName = getFoamName(longName);
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

        foamVtuData& vtuData = cachedVtu_(longName);

        vtkSmartPointer<vtkUnstructuredGrid> vtkgeom;
        if (vtuData.nPoints() && vtuData.pointMap().size())
        {
            if (meshState_ == polyMesh::UNCHANGED)
            {
                if (debug)
                {
                    Info<< "reuse " << longName << nl;
                }
                vtuData.reuse();
                continue;
            }
            else if (meshState_ == polyMesh::POINTS_MOVED)
            {
                if (debug)
                {
                    Info<< "move points " << longName << nl;
                }
                vtkgeom = vtuData.getCopy();
                vtkgeom->SetPoints
                (
                    movePoints(mesh, vtuData, vtuData.pointMap())
                );
            }
        }

        if (!vtkgeom)
        {
            fvMeshSubset subsetter(mesh);
            subsetter.setLargeCellSubset(zMesh[zoneId]);

            vtkgeom = volumeVTKSubsetMesh(subsetter, vtuData);
        }

        vtuData.set(vtkgeom);
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
        const word partName = getFoamName(longName);

        if (debug)
        {
            Info<< "Creating VTK mesh for cellSet=" << partName << endl;
        }

        foamVtuData& vtuData = cachedVtu_(longName);

        vtkSmartPointer<vtkUnstructuredGrid> vtkgeom;
        if (vtuData.nPoints() && vtuData.pointMap().size())
        {
            if (meshState_ == polyMesh::UNCHANGED)
            {
                if (debug)
                {
                    Info<< "reuse " << longName << nl;
                }
                vtuData.reuse();
                continue;
            }
            else if (meshState_ == polyMesh::POINTS_MOVED)
            {
                if (debug)
                {
                    Info<< "move points " << longName << nl;
                }
                vtkgeom = vtuData.getCopy();
                vtkgeom->SetPoints
                (
                    movePoints(mesh, vtuData, vtuData.pointMap())
                );
            }
        }

        if (!vtkgeom)
        {
            fvMeshSubset subsetter(mesh);
            subsetter.setLargeCellSubset(cellSet(mesh, partName));

            vtkgeom = volumeVTKSubsetMesh(subsetter, vtuData);
        }

        vtuData.set(vtkgeom);
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshFaceZones()
{
    const arrayRange& range = rangeFaceZones_;
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
        const word zoneName = getFoamName(longName);
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

        foamVtpData& vtpData = cachedVtp_(longName);

        vtkSmartPointer<vtkPolyData> vtkgeom;
        if (vtpData.nPoints())
        {
            if (meshState_ == polyMesh::UNCHANGED)
            {
                // Without movement is easy.
                if (debug)
                {
                    Info<<"reuse " << longName << nl;
                }
                vtpData.reuse();
                continue;
            }
            else if (meshState_ == polyMesh::POINTS_MOVED)
            {
                // Need point maps etc - not worth it at the moment
            }
        }

        if (!vtkgeom)
        {
            vtpData.clear();  // No additional ids, maps

            const primitiveFacePatch& pp = zMesh[zoneId]();
            vtkgeom = patchVTKMesh(pp);
        }

        vtpData.set(vtkgeom);
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshFaceSets()
{
    const arrayRange& range = rangeFaceSets_;
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
        const word partName = getFoamName(longName);

        if (debug)
        {
            Info<< "Creating VTK mesh for faceSet=" << partName << endl;
        }

        foamVtpData& vtpData = cachedVtp_(longName);

        vtkSmartPointer<vtkPolyData> vtkgeom;
        if (vtpData.nPoints())
        {
            if (meshState_ == polyMesh::UNCHANGED)
            {
                // Without movement is easy.
                if (debug)
                {
                    Info<<"reuse " << longName << nl;
                }
                vtpData.reuse();
                continue;
            }
            else if (meshState_ == polyMesh::POINTS_MOVED)
            {
                // Need point maps etc - not worth it at the moment
            }
        }

        if (!vtkgeom)
        {
            vtpData.clear();  // No other additional ids, maps

            // Misuse cellMap for face labels - sorted order for reliability
            vtpData.cellMap() = faceSet(mesh, partName).sortedToc();

            if (vtpData.cellMap().size())
            {
                uindirectPrimitivePatch pp
                (
                    UIndirectList<face>(mesh.faces(), vtpData.cellMap()),
                    mesh.points()
                );

                vtkgeom = patchVTKMesh(pp);
            }
        }

        if (vtkgeom)
        {
            vtpData.set(vtkgeom);
        }
        else
        {
            cachedVtp_.erase(longName);
        }
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertMeshPointZones()
{
    const arrayRange& range = rangePointZones_;
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
            const word zoneName = getFoamName(longName);
            const label zoneId = zMesh.findZoneID(zoneName);

            if (zoneId < 0)
            {
                continue;
            }

            foamVtpData& vtpData = cachedVtp_(longName);

            vtkSmartPointer<vtkPolyData> vtkgeom;
            if (vtpData.nPoints() && vtpData.pointMap().size())
            {
                if (meshState_ == polyMesh::UNCHANGED)
                {
                    if (debug)
                    {
                        Info<< "reusing " << longName << nl;
                    }
                    vtpData.reuse();
                    continue;
                }
                else if (meshState_ == polyMesh::POINTS_MOVED)
                {
                    if (debug)
                    {
                        Info<< "move points " << longName << nl;
                    }
                    vtkgeom = vtpData.getCopy();
                }
            }

            if (!vtkgeom)
            {
                // First time, or topo change
                vtkgeom = vtkSmartPointer<vtkPolyData>::New();
                vtpData.pointMap() = zMesh[zoneId];
            }

            const pointField& points = mesh.points();
            const labelUList& pointMap = vtpData.pointMap();

            auto vtkpoints = vtkSmartPointer<vtkPoints>::New();

            vtkpoints->SetNumberOfPoints(pointMap.size());
            forAll(pointMap, pointi)
            {
                vtkpoints->SetPoint(pointi, points[pointMap[pointi]].v_);
            }

            vtkgeom->SetPoints(vtkpoints);
            vtpData.set(vtkgeom);
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
    const arrayRange& range = rangePointSets_;
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
        const word partName = getFoamName(longName);

        foamVtpData& vtpData = cachedVtp_(longName);

        vtkSmartPointer<vtkPolyData> vtkgeom;
        if (vtpData.nPoints() && vtpData.pointMap().size())
        {
            if (meshState_ == polyMesh::UNCHANGED)
            {
                if (debug)
                {
                    Info<< "reusing " << longName << nl;
                }
                vtpData.reuse();
                continue;
            }
            else if (meshState_ == polyMesh::POINTS_MOVED)
            {
                if (debug)
                {
                    Info<< "move points " << longName << nl;
                }
                vtkgeom = vtpData.getCopy();
            }
        }

        if (!vtkgeom)
        {
            // First time, or topo change
            vtkgeom = vtkSmartPointer<vtkPolyData>::New();
            vtpData.pointMap() = pointSet(mesh, partName).sortedToc();
        }

        const pointField& points = mesh.points();
        const labelUList& pointMap = vtpData.pointMap();

        auto vtkpoints = vtkSmartPointer<vtkPoints>::New();

        vtkpoints->SetNumberOfPoints(pointMap.size());
        forAll(pointMap, pointi)
        {
            vtkpoints->SetPoint(pointi, points[pointMap[pointi]].v_);
        }

        vtkgeom->SetPoints(vtkpoints);
        vtpData.set(vtkgeom);
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << endl;
        printMemory();
    }
}


// ************************************************************************* //
