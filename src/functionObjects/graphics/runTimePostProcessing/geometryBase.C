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

#include "geometryBase.H"
#include "runTimePostProcessing.H"
#include "Constant.H"

#include "vtkActor.h"
#include "vtkProperty.h"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::functionObjects::runTimePostPro::geometryBase::renderModeType
>
Foam::functionObjects::runTimePostPro::geometryBase::renderModeTypeNames
{
    { renderModeType::rmFlat, "flat" },
    { renderModeType::rmGouraud, "gouraud" },
    { renderModeType::rmPhong, "phong" },
};


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::runTimePostPro::geometryBase::initialiseActor
(
    vtkActor* actor
) const
{
    actor->GetProperty()->SetSpecular(0);
    actor->GetProperty()->SetSpecularPower(20);

    switch (renderMode_)
    {
        case rmFlat:
        {
            actor->GetProperty()->SetInterpolationToFlat();
            break;
        }
        case rmGouraud:
        {
            actor->GetProperty()->SetInterpolationToGouraud();
            break;
        }
        case rmPhong:
        {
            actor->GetProperty()->SetInterpolationToPhong();
            break;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::geometryBase::geometryBase
(
    const runTimePostProcessing& parent,
    const dictionary& dict,
    const HashPtrTable<Function1<vector>, word>& colours
)
:
    parent_(parent),
    name_(dict.dictName()),
    visible_(readBool(dict.lookup("visible"))),
    renderMode_(rmGouraud),
    opacity_(nullptr),
    colours_(colours)
{
    if (dict.found("renderMode"))
    {
        renderMode_ = renderModeTypeNames.lookup("renderMode", dict);
    }

    if (dict.found("opacity"))
    {
        opacity_.reset(Function1<scalar>::New("opacity", dict).ptr());
    }
    else
    {
        opacity_.reset(new Function1Types::Constant<scalar>("opacity", 1.0));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::geometryBase::~geometryBase()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::functionObjects::runTimePostProcessing&
Foam::functionObjects::runTimePostPro::geometryBase::parent() const
{
    return parent_;
}


const Foam::word&
Foam::functionObjects::runTimePostPro::geometryBase::name() const
{
    return name_;
}


bool Foam::functionObjects::runTimePostPro::geometryBase::visible() const
{
    return visible_;
}


Foam::scalar Foam::functionObjects::runTimePostPro::geometryBase::opacity
(
    const scalar position
) const
{
    return opacity_->value(position);
}


const Foam::HashPtrTable<Foam::Function1<Foam::vector>, Foam::word>&
Foam::functionObjects::runTimePostPro::geometryBase::colours() const
{
    return colours_;
}


// ************************************************************************* //
