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

// VTK includes
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class PatchType>
vtkSmartPointer<vtkPolyData> Foam::vtkPVFoam::patchVTKMesh
(
    const string& name,
    const PatchType& p
)
{
    vtkSmartPointer<vtkPolyData> vtkmesh =
        vtkSmartPointer<vtkPolyData>::New();

    if (debug)
    {
        Info<< "<beg> patchVTKMesh - " << name << endl;
        printMemory();
    }

    // Convert OpenFOAM mesh vertices to VTK
    const Foam::pointField& points = p.localPoints();

    vtkSmartPointer<vtkPoints> vtkpoints =
        vtkSmartPointer<vtkPoints>::New();

    vtkpoints->SetNumberOfPoints(points.size());
    forAll(points, i)
    {
        vtkpoints->SetPoint(i, points[i].v_);
    }
    vtkmesh->SetPoints(vtkpoints);

    // Add faces as polygons
    const faceList& faces = p.localFaces();

    label nAlloc = faces.size();
    forAll(faces, facei)
    {
        nAlloc += faces[facei].size();
    }

    vtkSmartPointer<vtkCellArray> vtkcells =
        vtkSmartPointer<vtkCellArray>::New();

    vtkcells->Allocate(nAlloc);
    // If reusing memory, ensure insert always starts from 0
    vtkcells->Reset();

    forAll(faces, facei)
    {
        const face& f = faces[facei];
        vtkIdType nodeIds[f.size()];

        forAll(f, fp)
        {
            nodeIds[fp] = f[fp];
        }
        vtkcells->InsertNextCell(f.size(), nodeIds);
    }

    vtkmesh->SetPolys(vtkcells);

    if (debug)
    {
        Info<< "<end> patchVTKMesh - " << name << endl;
        printMemory();
    }

    return vtkmesh;
}


// ************************************************************************* //
