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
#include "pointData.H"
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
    const char* NamedEnum<pointData::representationType, 2>::names[] =
    {
        "sphere",
        "vector"
    };

    defineTypeNameAndDebug(pointData, 0);
    defineRunTimeSelectionTable(pointData, dictionary);
}

const Foam::NamedEnum<Foam::pointData::representationType, 2>
    Foam::pointData::representationTypeNames;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::pointData::addPoints
(
    const label frameI,
    vtkActor* actor,
    vtkPolyDataMapper* mapper,
    vtkPolyData* data
) const
{
    geometryBase::initialiseActor(actor);

    vector colour = pointColour_->value(frameI);
    actor->GetProperty()->SetColor(colour[0], colour[1], colour[2]);

    switch (representation_)
    {
        case rtSphere:
        case rtVector:
        {
            break;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointData::pointData
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
    maxGlyphLength_(readScalar(dict.lookup("maxGlyphLength"))),
    pointColour_(NULL)
{
    if (dict.found("pointColour"))
    {
        pointColour_.reset(DataEntry<vector>::New("pointColour", dict).ptr());
    }
    else
    {
        pointColour_.reset(colours["point"]->clone().ptr());
    }

    switch (representation_)
    {
        case rtSphere:
        {
            break;
        }
        case rtVector:
        {
            break;
        }
    }

}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::pointData> Foam::pointData::New
(
    const runTimePostProcessing& parent,
    const dictionary& dict,
    const HashPtrTable<DataEntry<vector>, word>& colours,
    const word& pointDataType
)
{
    if (debug)
    {
        Info<< "Selecting pointData " << pointDataType << endl;
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(pointDataType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "Foam::autoPtr<Foam::pointData> Foam::pointData::New"
            "("
                "const runTimePostProcessing&, "
                "const dictionary&, "
                "const HashPtrTable<DataEntry<vector>, word>&, "
                "const word&"
            ")"
        )   << "Unknown pointData type "
            << pointDataType << nl << nl
            << "Valid pointData types are:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<pointData>(cstrIter()(parent, dict, colours));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointData::~pointData()
{}


// ************************************************************************* //
