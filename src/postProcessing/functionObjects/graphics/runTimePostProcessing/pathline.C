/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  |
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
#include "pathline.H"
#include "runTimePostProcessing.H"

// VTK includes
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkSmartPointer.h"
#include "vtkTubeFilter.h"
#include "vtkLookupTable.h"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* NamedEnum<pathline::representationType, 4>::names[] =
    {
        "none",
        "line",
        "tube",
        "vector"
    };

    defineTypeNameAndDebug(pathline, 0);
    defineRunTimeSelectionTable(pathline, dictionary);
}

const Foam::NamedEnum<Foam::pathline::representationType, 4>
    Foam::pathline::representationTypeNames;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::pathline::addLines
(
    const label frameI,
    vtkActor* actor,
    vtkPolyData* data
) const
{
    geometryBase::initialiseActor(actor);

    vector colour = lineColour_->value(frameI);
    actor->GetProperty()->SetColor(colour[0], colour[1], colour[2]);

    vtkPolyDataMapper* mapper =
            vtkPolyDataMapper::SafeDownCast(actor->GetMapper());

    switch (representation_)
    {
        case rtNone:
        {
            actor->VisibilityOff();
            break;
        }
        case rtLine:
        {
            mapper->SetInputData(data);
            mapper->Update();
            break;

        }
        case rtTube:
        {
            vtkSmartPointer<vtkTubeFilter> tubes =
                vtkSmartPointer<vtkTubeFilter>::New();
            tubes->SetInputData(data);
            tubes->SetRadius(tubeRadius_);
            tubes->SetNumberOfSides(20);
            tubes->CappingOn();
            tubes->Update();

            mapper->SetInputConnection(tubes->GetOutputPort());
            mapper->Update();

            break;

        }
        case rtVector:
        {
            break;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pathline::pathline
(
    const runTimePostProcessing& parent,
    const dictionary& dict,
    const HashPtrTable<DataEntry<vector>, word>& colours
)
:
    geometryBase(parent, dict, colours),
    representation_
    (
        representationTypeNames.read(dict.lookup("representation"))
    ),
    tubeRadius_(0.0),
    lineColour_(NULL)
{
    if (dict.found("lineColour"))
    {
        lineColour_.reset(DataEntry<vector>::New("lineColour", dict).ptr());
    }
    else
    {
        lineColour_.reset(colours["line"]->clone().ptr());
    }

    switch (representation_)
    {
        case rtNone:
        {
            break;
        }
        case rtLine:
        {
            break;
        }
        case rtTube:
        {
            dict.lookup("tubeRadius") >> tubeRadius_;
            break;
        }
        case rtVector:
        {
            break;
        }
    }

}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::pathline> Foam::pathline::New
(
    const runTimePostProcessing& parent,
    const dictionary& dict,
    const HashPtrTable<DataEntry<vector>, word>& colours,
    const word& pathlineType
)
{
    if (debug)
    {
        Info<< "Selecting pathline " << pathlineType << endl;
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(pathlineType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "Foam::autoPtr<Foam::pathline> Foam::pathline::New"
            "("
                "const runTimePostProcessing&, "
                "const dictionary&, "
                "const HashPtrTable<DataEntry<vector>, word>&, "
                "const word&"
            ")"
        )   << "Unknown pathline type "
            << pathlineType << nl << nl
            << "Valid pathline types are:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<pathline>(cstrIter()(parent, dict, colours));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pathline::~pathline()
{}


// ************************************************************************* //
