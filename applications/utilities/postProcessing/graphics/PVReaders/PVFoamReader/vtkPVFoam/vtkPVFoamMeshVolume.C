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
#include "foamVtuSizing.H"

// VTK includes
#include "vtkUnstructuredGrid.h"
#include "foamVtkAdaptors.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vtkSmartPointer<vtkUnstructuredGrid> Foam::vtkPVFoam::volumeVTKMesh
(
    const fvMesh& mesh,
    foamVtuData& vtuData
)
{
    if (debug)
    {
        Info<< "<beg> " << FUNCTION_NAME << endl;
        printMemory();
    }

    foamVtuSizing sizing(mesh, !reader_->GetUseVTKPolyhedron());

    vtkSmartPointer<vtkUnstructuredGrid> vtkmesh =
        vtkSmartPointer<vtkUnstructuredGrid>::New();

    auto cellTypes =
        nonNullSmartPointer<vtkUnsignedCharArray>(vtkmesh->GetCellTypesArray());

    auto cells =
        nonNullSmartPointer<vtkCellArray>(vtkmesh->GetCells());

    auto faces =
        nonNullSmartPointer<vtkIdTypeArray>(vtkmesh->GetFaces());

    auto cellLocations =
        nonNullSmartPointer<vtkIdTypeArray>(vtkmesh->GetCellLocationsArray());

    auto faceLocations =
        nonNullSmartPointer<vtkIdTypeArray>(vtkmesh->GetFaceLocations());

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
            sizing.sizeInternal(foamVtuSizing::slotType::CELLS)
        );

    UList<vtkIdType> cellLocationsUL =
        vtkUList
        (
            cellLocations,
            sizing.sizeInternal(foamVtuSizing::slotType::CELLS_OFFSETS)
        );

    UList<vtkIdType> facesUL =
        vtkUList
        (
            faces,
            sizing.sizeInternal(foamVtuSizing::slotType::FACES)
        );

    UList<vtkIdType> faceLocationsUL =
        vtkUList
        (
            faceLocations,
            sizing.sizeInternal(foamVtuSizing::slotType::FACES_OFFSETS)
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


    // Convert OpenFOAM mesh vertices to VTK
    // - can only do this *after* populating the decompInfo with cell-ids
    //   for any additional points (ie, mesh cell-centres)
    vtkSmartPointer<vtkPoints> vtkpoints = vtkmesh->GetPoints();
    if (!vtkpoints)
    {
        // No points previously, add an entry
        vtkpoints = vtkSmartPointer<vtkPoints>::New();
        vtkmesh->SetPoints(vtkpoints);
    }

    vtkpoints->SetNumberOfPoints(sizing.nFieldPoints());

    // Normal points
    const pointField& points = mesh.points();
    const labelList& addPoints = vtuData.additionalIds();

    label pointi = 0;
    forAll(points, i)
    {
        vtkpoints->SetPoint(pointi++, points[i].v_);
    }
    forAll(addPoints, i)
    {
        vtkpoints->SetPoint(pointi++, mesh.C()[addPoints[i]].v_);
    }

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
        Info<< "<end> " << FUNCTION_NAME << endl;
        printMemory();
    }

    return vtkmesh;
}


// ************************************************************************* //
