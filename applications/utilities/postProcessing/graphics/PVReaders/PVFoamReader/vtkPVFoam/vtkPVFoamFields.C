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
#include "Cloud.H"
#include "IOobjectList.H"
#include "vtkPVFoamReader.h"

// VTK includes
#include "vtkDataArraySelection.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"

// Templates (only needed here)
#include "vtkPVFoamFieldTemplates.C"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPVFoam::pruneObjectList
(
    IOobjectList& objects,
    const hashedWordList& retain
)
{
    if (retain.empty())
    {
        objects.clear();
    }

    // only retain specified fields
    forAllIters(objects, iter)
    {
        if (!retain.found(iter()->name()))
        {
            objects.erase(iter);
        }
    }
}


void Foam::vtkPVFoam::convertVolFields()
{
    const fvMesh& mesh = *meshPtr_;

    const bool interpFields = reader_->GetInterpolateVolFields();
    hashedWordList selectedFields = getSelected
    (
        reader_->GetVolFieldSelection()
    );

    if (selectedFields.empty())
    {
        return;
    }

    // Get objects (fields) for this time - only keep selected fields
    // the region name is already in the mesh db
    IOobjectList objects(mesh, dbPtr_().timeName());
    pruneObjectList(objects, selectedFields);

    if (objects.empty())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> convert volume fields" << endl;
        forAllConstIters(objects, iter)
        {
            Info<< "  " << iter()->name()
                << " == " << iter()->objectPath() << endl;
        }
        printMemory();
    }


    PtrList<patchInterpolator> interpLst;

    if (interpFields)
    {
        interpLst.setSize(mesh.boundaryMesh().size());

        forAll(interpLst, i)
        {
            interpLst.set
            (
                i,
                new PrimitivePatchInterpolation<primitivePatch>
                (
                    mesh.boundaryMesh()[i]
                )
            );
        }
    }

    convertVolFields<scalar>(mesh, interpLst, objects);
    convertVolFields<vector>(mesh, interpLst, objects);
    convertVolFields<sphericalTensor>(mesh, interpLst, objects);
    convertVolFields<symmTensor>(mesh, interpLst, objects);
    convertVolFields<tensor>(mesh, interpLst, objects);

    convertDimFields<scalar>(mesh, interpLst, objects);
    convertDimFields<vector>(mesh, interpLst, objects);
    convertDimFields<sphericalTensor>(mesh, interpLst, objects);
    convertDimFields<symmTensor>(mesh, interpLst, objects);
    convertDimFields<tensor>(mesh, interpLst, objects);

    if (debug)
    {
        Info<< "<end> convert volume fields" << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertPointFields()
{
    const fvMesh& mesh = *meshPtr_;

    hashedWordList selectedFields = getSelected
    (
        reader_->GetPointFieldSelection()
    );

    if (selectedFields.empty())
    {
        if (debug)
        {
            Info<< "no point fields selected" << endl;
        }
        return;
    }

    // Get objects (fields) for this time - only keep selected fields
    // the region name is already in the mesh db
    IOobjectList objects(mesh, dbPtr_().timeName());
    pruneObjectList(objects, selectedFields);

    if (objects.empty())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> convert volume -> point fields" << endl;
        forAllConstIters(objects, iter)
        {
            Info<< "  " << iter()->name()
                << " == " << iter()->objectPath() << endl;
        }
        printMemory();
    }

    // Construct interpolation on the raw mesh
    const pointMesh& pMesh = pointMesh::New(mesh);

    convertPointFields<scalar>(pMesh, objects);
    convertPointFields<vector>(pMesh, objects);
    convertPointFields<sphericalTensor>(pMesh, objects);
    convertPointFields<symmTensor>(pMesh, objects);
    convertPointFields<tensor>(pMesh, objects);

    if (debug)
    {
        Info<< "<end> convert volume -> point fields" << endl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertLagrangianFields()
{
    arrayRange& range = rangeLagrangian_;
    const fvMesh& mesh = *meshPtr_;

    hashedWordList selectedFields = getSelected
    (
        reader_->GetLagrangianFieldSelection()
    );

    if (selectedFields.empty())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> convert Lagrangian fields" << endl;
        printMemory();
    }

    for (auto partId : range)
    {
        if (!selectedPartIds_.found(partId))
        {
            continue;
        }

        const word cloudName = getPartName(partId);
        auto iter = cachedVtp_.find(selectedPartIds_[partId]);

        if (!iter.found() || !(*iter).vtkmesh)
        {
            continue;
        }
        auto vtkmesh = (*iter).vtkmesh;

        // Get the Lagrangian fields for this time and this cloud
        // but only keep selected fields
        // the region name is already in the mesh db
        IOobjectList objects
        (
            mesh,
            dbPtr_().timeName(),
            cloud::prefix/cloudName
        );
        pruneObjectList(objects, selectedFields);

        if (objects.empty())
        {
            continue;
        }

        if (debug)
        {
            Info<< "converting OpenFOAM lagrangian fields" << endl;
            forAllConstIters(objects, iter)
            {
                Info<< "  " << iter()->name()
                    << " == " << iter()->objectPath() << endl;
            }
        }

        convertLagrangianFields<label>(objects, vtkmesh);
        convertLagrangianFields<scalar>(objects, vtkmesh);
        convertLagrangianFields<vector>(objects, vtkmesh);
        convertLagrangianFields<sphericalTensor>(objects, vtkmesh);
        convertLagrangianFields<symmTensor>(objects, vtkmesh);
        convertLagrangianFields<tensor>(objects, vtkmesh);
    }

    if (debug)
    {
        Info<< "<end> convert Lagrangian fields" << endl;
        printMemory();
    }
}


// ************************************************************************* //
