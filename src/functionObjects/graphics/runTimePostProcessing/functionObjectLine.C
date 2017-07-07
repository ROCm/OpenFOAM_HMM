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
#include "functionObjectLine.H"
#include "runTimePostProcessing.H"
#include "addToRunTimeSelectionTable.H"

// VTK includes
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataReader.h"
#include "vtkProperty.h"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{
    defineTypeNameAndDebug(functionObjectLine, 0);
    addToRunTimeSelectionTable(pathline, functionObjectLine, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::functionObjectLine::functionObjectLine
(
    const runTimePostProcessing& parent,
    const dictionary& dict,
    const HashPtrTable<Function1<vector>, word>& colours
)
:
    pathline(parent, dict, colours),
    functionObjectBase(parent, dict, colours),
    actor_()
{
    actor_ = vtkSmartPointer<vtkActor>::New();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::functionObjectLine::~functionObjectLine()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::functionObjects::runTimePostPro::functionObjectLine::
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
            << ". Line will not be processed"
            << endl;
        return;
    }

    if (fName.ext() == "vtk")
    {
        auto lines = vtkSmartPointer<vtkPolyDataReader>::New();
        lines->SetFileName(fName.c_str());
        lines->Update();

        auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        setField(position, fieldName_, mapper, renderer, lines->GetOutput());

        actor_->SetMapper(mapper);

        addLines(position, actor_, lines->GetOutput());

        renderer->AddActor(actor_);
    }
}


void Foam::functionObjects::runTimePostPro::functionObjectLine::updateActors
(
    const scalar position
)
{
    actor_->GetProperty()->SetLineWidth(2);
    actor_->GetProperty()->SetOpacity(opacity(position));
}

bool Foam::functionObjects::runTimePostPro::functionObjectLine::clear()
{
    if (functionObjectBase::clear())
    {
        return removeFile("file", fieldName_);
    }

    return false;
}


// ************************************************************************* //
