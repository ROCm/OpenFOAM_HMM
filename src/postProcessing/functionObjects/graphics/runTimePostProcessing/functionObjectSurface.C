/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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

// OpenFOAM includes
#include "functionObjectSurface.H"
#include "runTimePostProcessing.H"
#include "addToRunTimeSelectionTable.H"

// VTK includes
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataReader.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(functionObjectSurface, 0);
    addToRunTimeSelectionTable(surface, functionObjectSurface, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjectSurface::functionObjectSurface
(
    const runTimePostProcessing& parent,
    const dictionary& dict,
    const HashPtrTable<Function1<vector>, word>& colours
)
:
    geometrySurface(parent, dict, colours, List<fileName>()),
    fieldVisualisationBase(parent, dict, colours),
    functionObject_("")
{
    if (visible_)
    {
        dict.lookup("functionObject") >> functionObject_;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjectSurface::~functionObjectSurface()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::functionObjectSurface::addGeometryToScene
(
    const scalar position,
    vtkRenderer* renderer
)
{
    if (!visible_)
    {
        return;
    }

    const dictionary dict =
        geometryBase::parent_.getObjectProperty
        (
            functionObject_,
            fieldName_,
            dictionary::null
        );

    fileName fName;
    if (!dict.readIfPresent("file", fName))
    {
        WarningInFunction
            << "Unable to find function object " << functionObject_
            << " output for field " << fieldName_
            << ". Surface will not be processed"
            << endl;
        return;
    }


    if (representation_ == rtGlyph)
    {
        vtkSmartPointer<vtkPolyDataReader> surf =
            vtkSmartPointer<vtkPolyDataReader>::New();
        surf->SetFileName(fName.c_str());
        surf->Update();

        addGlyphs
        (
            position,
            fieldName_,
            fieldName_,
            maxGlyphLength_,
            surf->GetOutput(),
            surfaceActor_,
            renderer
        );
    }
    else
    {
        if (fName.ext() == "vtk")
        {
            vtkSmartPointer<vtkPolyDataReader> surf =
                vtkSmartPointer<vtkPolyDataReader>::New();
            surf->SetFileName(fName.c_str());
            surf->Update();

            addFeatureEdges(renderer, surf->GetOutput());

            vtkSmartPointer<vtkPolyDataMapper> mapper =
                vtkSmartPointer<vtkPolyDataMapper>::New();
            mapper->SetInputConnection(surf->GetOutputPort());

            setField(position, fieldName_, mapper, renderer);

            surfaceActor_->SetMapper(mapper);

            setRepresentation(surfaceActor_);

            renderer->AddActor(surfaceActor_);
        }
        else
        {
            WarningInFunction
                << "Only VTK file types are supported"
                << endl;
        }
    }
}


// ************************************************************************* //
