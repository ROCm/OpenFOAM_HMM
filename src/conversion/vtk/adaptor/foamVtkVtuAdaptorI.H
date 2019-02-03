/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2019 OpenCFD Ltd.
     \\/     M anipulation  |
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline vtkSmartPointer<vtkPoints>
Foam::vtk::vtuAdaptor::points
(
    const fvMesh& mesh
) const
{
    // Convert OpenFOAM mesh vertices to VTK
    auto vtkpoints = vtkSmartPointer<vtkPoints>::New();

    // Normal points
    const pointField& pts = mesh.points();

    // Additional cell centres
    const labelUList& addPoints = this->additionalIds();

    vtkpoints->SetNumberOfPoints(pts.size() + addPoints.size());

    // Normal points
    vtkIdType pointId = 0;
    for (const point& p : pts)
    {
        vtkpoints->SetPoint(pointId++, p.v_);
    }

    // Cell centres
    for (const label meshCelli : addPoints)
    {
        vtkpoints->SetPoint(pointId++, mesh.cellCentres()[meshCelli].v_);
    }

    return vtkpoints;
}


inline vtkSmartPointer<vtkPoints>
Foam::vtk::vtuAdaptor::points
(
    const fvMesh& mesh,
    const labelUList& pointMap
) const
{
    // Convert OpenFOAM mesh vertices to VTK
    auto vtkpoints = vtkSmartPointer<vtkPoints>::New();

    // Normal points
    const pointField& pts = mesh.points();

    // Additional cell centres
    const labelUList& addPoints = this->additionalIds();

    vtkpoints->SetNumberOfPoints(pointMap.size() + addPoints.size());

    // Normal points
    vtkIdType pointId = 0;
    for (const label meshPointi : pointMap)
    {
        vtkpoints->SetPoint(pointId++, pts[meshPointi].v_);
    }

    // Cell centres
    for (const label meshCelli : addPoints)
    {
        vtkpoints->SetPoint(pointId++, mesh.cellCentres()[meshCelli].v_);
    }

    return vtkpoints;
}


inline vtkSmartPointer<vtkUnstructuredGrid>
Foam::vtk::vtuAdaptor::internal
(
    const fvMesh& mesh,
    const bool decompPoly
)
{
    vtk::vtuSizing sizing(mesh, decompPoly);

    auto cellTypes = vtkSmartPointer<vtkUnsignedCharArray>::New();

    auto cells = vtkSmartPointer<vtkCellArray>::New();
    auto faces = vtkSmartPointer<vtkIdTypeArray>::New();

    auto cellLocations = vtkSmartPointer<vtkIdTypeArray>::New();
    auto faceLocations = vtkSmartPointer<vtkIdTypeArray>::New();

    UList<uint8_t> cellTypesUL =
        vtk::Tools::asUList(cellTypes, sizing.nFieldCells());

    UList<vtkIdType> cellsUL =
        vtk::Tools::asUList
        (
            cells,
            sizing.nFieldCells(),
            sizing.sizeInternal(vtk::vtuSizing::slotType::CELLS)
        );

    UList<vtkIdType> cellLocationsUL =
        vtk::Tools::asUList
        (
            cellLocations,
            sizing.sizeInternal(vtk::vtuSizing::slotType::CELLS_OFFSETS)
        );

    UList<vtkIdType> facesUL =
        vtk::Tools::asUList
        (
            faces,
            sizing.sizeInternal(vtk::vtuSizing::slotType::FACES)
        );

    UList<vtkIdType> faceLocationsUL =
        vtk::Tools::asUList
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
        static_cast<foamVtkMeshMaps&>(*this)
    );

    auto vtkmesh = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // Convert OpenFOAM mesh vertices to VTK
    // - can only do this *after* populating the decompInfo with cell-ids
    //   for any additional points (ie, mesh cell-centres)
    vtkmesh->SetPoints(this->points(mesh));

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

    return vtkmesh;
}


// ************************************************************************* //
