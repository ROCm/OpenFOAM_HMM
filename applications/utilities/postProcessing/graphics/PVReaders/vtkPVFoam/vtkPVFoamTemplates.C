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
#include "polyPatch.H"
#include "primitivePatch.H"
#include "foamVtkAdaptors.H"

// VTK includes
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class PatchType>
vtkSmartPointer<vtkPoints> Foam::vtkPVFoam::movePatchPoints
(
    const PatchType& p
)
{
    // Convert OpenFOAM mesh vertices to VTK
    const pointField& points = p.localPoints();

    auto vtkpoints = vtkSmartPointer<vtkPoints>::New();

    vtkpoints->SetNumberOfPoints(points.size());
    forAll(points, i)
    {
        vtkpoints->SetPoint(i, points[i].v_);
    }

    return vtkpoints;
}


template<class PatchType>
vtkSmartPointer<vtkCellArray> Foam::vtkPVFoam::patchFacesVTKCells
(
    const PatchType& p
)
{
    // Faces as polygons
    const faceList& faces = p.localFaces();

    label nAlloc = faces.size();
    for (const face& f : faces)
    {
        nAlloc += f.size();
    }

    auto cells = vtkSmartPointer<vtkCellArray>::New();

    UList<vtkIdType> cellsUL =
        vtkUList
        (
            cells,
            faces.size(),
            nAlloc
        );

    // Cell connectivity for polygons
    // [size, verts..., size, verts... ]
    label idx = 0;
    for (const face& f : faces)
    {
        cellsUL[idx++] = f.size();

        for (const label verti : f)
        {
            cellsUL[idx++] = verti;
        }
    }

    return cells;
}


template<class PatchType>
vtkSmartPointer<vtkPolyData> Foam::vtkPVFoam::patchVTKMesh
(
    const PatchType& p
)
{
    auto vtkmesh = vtkSmartPointer<vtkPolyData>::New();

    vtkmesh->SetPoints(movePatchPoints(p));
    vtkmesh->SetPolys(patchFacesVTKCells(p));

    return vtkmesh;
}


// ************************************************************************* //
