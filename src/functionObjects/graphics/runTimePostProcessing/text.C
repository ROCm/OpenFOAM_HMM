/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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
#include "text.H"
#include "fvMesh.H"
#include "runTimePostProcessing.H"

// VTK includes
#include "vtkRenderer.h"
#include "vtkSmartPointer.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::text::text
(
    const runTimePostProcessing& parent,
    const dictionary& dict,
    const HashPtrTable<Function1<vector>, word>& colours
)
:
    geometryBase(parent, dict, colours),
    string_(dict.lookup("string")),
    position_(dict.lookup("position")),
    size_(readScalar(dict.lookup("size"))),
    colour_(nullptr),
    bold_(readBool(dict.lookup("bold"))),
    timeStamp_(dict.lookupOrDefault<bool>("timeStamp", false))
{
    if (dict.found("colour"))
    {
        colour_.reset(Function1<vector>::New("colour", dict).ptr());
    }
    else
    {
        colour_.reset(colours["text"]->clone().ptr());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::text::~text()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::functionObjects::runTimePostPro::text::addGeometryToScene
(
    const scalar position,
    vtkRenderer* renderer
)
{
    if (!visible_)
    {
        return;
    }

    auto actor = vtkSmartPointer<vtkTextActor>::New();

    // Concatenate string with timeStamp if true
    string textAndTime = string_;
    if (timeStamp_)
    {
        textAndTime =
            textAndTime + " " + geometryBase::parent_.mesh().time().timeName();
    }
    actor->SetInput(textAndTime.c_str());
    actor->GetTextProperty()->SetFontFamilyToArial();
    actor->GetTextProperty()->SetFontSize(size_);
    actor->GetTextProperty()->SetJustificationToLeft();
    actor->GetTextProperty()->SetVerticalJustificationToBottom();
    actor->GetTextProperty()->SetBold(bold_);

    const vector colour = colour_->value(position);
    actor->GetTextProperty()->SetColor(colour[0], colour[1], colour[2]);
    actor->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
    actor->GetPositionCoordinate()->SetValue
    (
        position_.first(),
        position_.second()
    );

    renderer->AddActor2D(actor);
}


void Foam::functionObjects::runTimePostPro::text::updateActors
(
    const scalar position
)
{
    // Do nothing - all handled by addGeometryToScene
}


bool Foam::functionObjects::runTimePostPro::text::clear()
{
    return true;
}


// ************************************************************************* //
