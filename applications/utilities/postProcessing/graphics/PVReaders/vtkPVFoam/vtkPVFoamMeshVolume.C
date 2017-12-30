/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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
#include "vtkPVFoamReader.h"

// OpenFOAM includes
#include "fvMesh.H"
#include "fvMeshSubset.H"
#include "foamVtkAdaptors.H"
#include "foamVtuSizing.H"

// VTK includes
#include "vtkUnstructuredGrid.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vtkSmartPointer<vtkPoints> Foam::vtkPVFoam::movePoints
(
    const fvMesh& mesh,
    const foamVtuData& vtuData
)
{
    // Convert OpenFOAM mesh vertices to VTK
    auto vtkpoints = vtkSmartPointer<vtkPoints>::New();

    // Normal points
    const pointField& points = mesh.points();

    // Additional cell centres
    const labelList& addPoints = vtuData.additionalIds();

    vtkpoints->SetNumberOfPoints(points.size() + addPoints.size());

    // Normal points
    label pointi = 0;
    forAll(points, i)
    {
        vtkpoints->SetPoint(pointi++, points[i].v_);
    }

    // Cell centres
    forAll(addPoints, i)
    {
        vtkpoints->SetPoint(pointi++, mesh.C()[addPoints[i]].v_);
    }

    return vtkpoints;
}


vtkSmartPointer<vtkPoints> Foam::vtkPVFoam::movePoints
(
    const fvMesh& mesh,
    const foamVtuData& vtuData,
    const labelUList& pointMap
)
{
    // Convert OpenFOAM mesh vertices to VTK
    auto vtkpoints = vtkSmartPointer<vtkPoints>::New();

    // Normal points
    const pointField& points = mesh.points();

    // Additional cell centres
    const labelList& addPoints = vtuData.additionalIds();

    vtkpoints->SetNumberOfPoints(pointMap.size() + addPoints.size());

    // Normal points
    label pointi = 0;
    forAll(pointMap, i)
    {
        vtkpoints->SetPoint(pointi++, points[pointMap[i]].v_);
    }

    // Cell centres
    forAll(addPoints, i)
    {
        vtkpoints->SetPoint(pointi++, mesh.C()[addPoints[i]].v_);
    }

    return vtkpoints;
}


vtkSmartPointer<vtkUnstructuredGrid> Foam::vtkPVFoam::volumeVTKMesh
(
    const fvMesh& mesh,
    foamVtuData& vtuData,
    const bool decompPoly
)
{
    if (debug)
    {
        Info<< "<beg> " << FUNCTION_NAME << nl;
        printMemory();
    }

    vtk::vtuSizing sizing(mesh, decompPoly);

    auto cellTypes = vtkSmartPointer<vtkUnsignedCharArray>::New();

    auto cells = vtkSmartPointer<vtkCellArray>::New();
    auto faces = vtkSmartPointer<vtkIdTypeArray>::New();

    auto cellLocations = vtkSmartPointer<vtkIdTypeArray>::New();
    auto faceLocations = vtkSmartPointer<vtkIdTypeArray>::New();

    UList<uint8_t> cellTypesUL =
        vtkUList
        (
            cellTypes,
            sizing.nFieldCells()
        );

    UList<vtkIdType> cellsUL =
        vtkUList
        (
            cells,
            sizing.nFieldCells(),
            sizing.sizeInternal(vtk::vtuSizing::slotType::CELLS)
        );

    UList<vtkIdType> cellLocationsUL =
        vtkUList
        (
            cellLocations,
            sizing.sizeInternal(vtk::vtuSizing::slotType::CELLS_OFFSETS)
        );

    UList<vtkIdType> facesUL =
        vtkUList
        (
            faces,
            sizing.sizeInternal(vtk::vtuSizing::slotType::FACES)
        );

    UList<vtkIdType> faceLocationsUL =
        vtkUList
        (
            faceLocations,
            sizing.sizeInternal(vtk::vtuSizing::slotType::FACES_OFFSETS)
        );


    sizing.populateInternal
    (
        mesh,
        cellTypesUL,
        cellsUL,
        cellLocationsUL,
        facesUL,
        faceLocationsUL,
        static_cast<foamVtkMeshMaps&>(vtuData)
    );

    auto vtkmesh = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // Convert OpenFOAM mesh vertices to VTK
    // - can only do this *after* populating the decompInfo with cell-ids
    //   for any additional points (ie, mesh cell-centres)
    vtkmesh->SetPoints(movePoints(mesh, vtuData));

    if (facesUL.size())
    {
        vtkmesh->SetCells
        (
            cellTypes,
            cellLocations,
            cells,
            faceLocations,
            faces
        );
    }
    else
    {
        vtkmesh->SetCells
        (
            cellTypes,
            cellLocations,
            cells,
            nullptr,
            nullptr
        );
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << nl;
        printMemory();
    }

    return vtkmesh;
}


vtkSmartPointer<vtkUnstructuredGrid> Foam::vtkPVFoam::volumeVTKSubsetMesh
(
    const fvMeshSubset& subsetter,
    foamVtuData& vtuData,
    const bool decompPoly
)
{
    vtkSmartPointer<vtkUnstructuredGrid> vtkmesh = volumeVTKMesh
    (
        subsetter.subMesh(),
        vtuData,
        decompPoly
    );

    // Convert cellMap, addPointCellLabels to global cell ids
    vtuData.renumberCells(subsetter.cellMap());

    // Copy pointMap as well, otherwise pointFields fail
    vtuData.pointMap() = subsetter.pointMap();

    return vtkmesh;
}


vtkSmartPointer<vtkUnstructuredGrid> Foam::vtkPVFoam::volumeVTKMesh
(
    const fvMesh& mesh,
    foamVtuData& vtuData
) const
{
    return volumeVTKMesh(mesh, vtuData, this->decomposePoly_);
}


vtkSmartPointer<vtkUnstructuredGrid> Foam::vtkPVFoam::volumeVTKSubsetMesh
(
    const fvMeshSubset& subsetter,
    foamVtuData& vtuData
) const
{
    return volumeVTKSubsetMesh(subsetter, vtuData, this->decomposePoly_);
}


// ************************************************************************* //
