/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description

\*---------------------------------------------------------------------------*/

#include "vtkPV3Foam.H"

// Foam includes
#include "polyPatch.H"
#include "vtkPV3FoamInsertNextPoint.H"

// VTK includes
#include "vtkCellArray.h"
#include "vtkUnstructuredGrid.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPV3Foam::addPatchMesh
(
    const polyPatch& p,
    vtkUnstructuredGrid *vtkPatch
)
{
    if (debug)
    {
        Info<< "Adding patch: " << p.name() << endl;
    }

    SetName(vtkPatch, p.name().c_str());

    if (debug)
    {
        Info<< "converting points" << endl;
    }

    const Foam::pointField& points = p.localPoints();

    // Convert Foam mesh vertices to VTK
    vtkPoints *vtkpoints = vtkPoints::New();
    vtkpoints->Allocate(points.size());

    forAll(points, i)
    {
        vtkPV3FoamInsertNextPoint(vtkpoints, points[i]);
    }

    if (debug)
    {
        Info<< "converting faces" << endl;
    }

    const faceList& faces = p.localFaces();

    vtkPatch->Allocate(faces.size());

    forAll(faces, facei)
    {
        const face& f = faces[facei];
        vtkIdType nodeIds[f.size()];

        if (f.size() == 3)
        {
            for (int j=0; j<3; j++)
            {
                nodeIds[j] = f[j];
            }
            vtkPatch->InsertNextCell
            (
                VTK_TRIANGLE,
                3,
                nodeIds
            );
        }
        else if (f.size() == 4)
        {
            for (int j=0; j<4; j++)
            {
                nodeIds[j] = f[j];
            }
            vtkPatch->InsertNextCell
            (
                VTK_QUAD,
                4,
                nodeIds
            );
        }
        else
        {
            for (int j=0; j<f.size(); j++)
            {
                nodeIds[j] = f[j];
            }
            vtkPatch->InsertNextCell
            (
                VTK_POLYGON,
                f.size(),
                nodeIds
            );
        }
    }


    vtkPatch->SetPoints(vtkpoints);
    vtkpoints->Delete();
}


// ************************************************************************* //
