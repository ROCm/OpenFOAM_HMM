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
namespace functionObjects
{
namespace runTimePostPro
{
    defineTypeNameAndDebug(functionObjectSurface, 0);
    addToRunTimeSelectionTable(surface, functionObjectSurface, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::functionObjectSurface::
functionObjectSurface
(
    const runTimePostProcessing& parent,
    const dictionary& dict,
    const HashPtrTable<Function1<vector>, word>& colours
)
:
    geometrySurface(parent, dict, colours, List<fileName>()),
    functionObjectBase(parent, dict, colours)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::functionObjectSurface::
~functionObjectSurface()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::functionObjects::runTimePostPro::functionObjectSurface::
addGeometryToScene
(
    const scalar position,
    vtkRenderer* renderer
)
{
    if (!visible_)
    {
        return;
    }

    fileName fName = getFileName("file", fieldName_);
    if (fName.empty())
    {
        WarningInFunction
            << "Unable to read file name from function object "
            << functionObjectName_ << " for field " << fieldName_
            << ". Surface will not be processed"
            << endl;
        return;
    }


    if (representation_ == rtGlyph)
    {
        auto surf = vtkSmartPointer<vtkPolyDataReader>::New();
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
            auto surf = vtkSmartPointer<vtkPolyDataReader>::New();
            surf->SetFileName(fName.c_str());
            surf->Update();

            addFeatureEdges(renderer, surf->GetOutput());

            auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            mapper->SetInputConnection(surf->GetOutputPort());

            setField(position, fieldName_, mapper, renderer, surf->GetOutput());

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


bool Foam::functionObjects::runTimePostPro::functionObjectSurface::clear()
{
    if (functionObjectBase::clear())
    {
        return removeFile("file", fieldName_);
    }

    return false;
}


// ************************************************************************* //
