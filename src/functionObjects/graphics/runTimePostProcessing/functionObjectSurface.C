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
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkSmartPointer.h"

// VTK Readers
#include "vtkPolyDataReader.h"
#include "vtkXMLPolyDataReader.h"

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


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

static vtkSmartPointer<vtkPolyData> getPolyDataFile(const Foam::fileName& fName)
{
    // Not extremely elegant...
    vtkSmartPointer<vtkPolyData> dataset;

    if (fName.ext() == "vtk")
    {
        auto reader = vtkSmartPointer<vtkPolyDataReader>::New();

        reader->SetFileName(fName.c_str());
        reader->Update();
        dataset = reader->GetOutput();

        return dataset;
    }

    if (fName.ext() == "vtp")
    {
        auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();

        reader->SetFileName(fName.c_str());
        reader->Update();
        dataset = reader->GetOutput();

        return dataset;
    }

    return dataset;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::functionObjectSurface::
functionObjectSurface
(
    const runTimePostProcessing& parent,
    const dictionary& dict,
    const HashPtrTable<Function1<vector>>& colours
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


    auto polyData = getPolyDataFile(fName);

    if (!polyData || polyData->GetNumberOfPoints() == 0)
    {
        WarningInFunction
            << "Could not read "<< fName << nl
            << "Only VTK (.vtp, .vtk) files are supported"
            << endl;
        return;
    }

    if (representation_ == rtGlyph)
    {
        addGlyphs
        (
            position,
            fieldName_,
            fieldName_,
            maxGlyphLength_,
            polyData,
            surfaceActor_,
            renderer
        );
    }
    else
    {
        addFeatureEdges(renderer, polyData);

        auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputData(polyData);

        setField(position, fieldName_, mapper, renderer, polyData);

        surfaceActor_->SetMapper(mapper);

        setRepresentation(surfaceActor_);

        renderer->AddActor(surfaceActor_);
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
