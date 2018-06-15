/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2018 OpenCFD Ltd.
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
#include "functionObjectCloud.H"
#include "fvMesh.H"
#include "runTimePostProcessing.H"
#include "addToRunTimeSelectionTable.H"
#include "stringOps.H"

// VTK includes
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkPolyDataReader.h"
#include "vtkProperty.h"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{
    defineTypeNameAndDebug(functionObjectCloud, 0);
    addToRunTimeSelectionTable(pointData, functionObjectCloud, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::functionObjectCloud::functionObjectCloud
(
    const runTimePostProcessing& parent,
    const dictionary& dict,
    const HashPtrTable<Function1<vector>, word>& colours
)
:
    pointData(parent, dict, colours),
    functionObjectBase(parent, dict, colours),
    cloudName_(dict.lookup("cloud")),
    inputFileName_(),
    colourFieldName_(dict.lookup("colourField")),
    actor_()
{
    actor_ = vtkSmartPointer<vtkActor>::New();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::functionObjectCloud::
~functionObjectCloud()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::functionObjects::runTimePostPro::functionObjectCloud::
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

    inputFileName_.clear();

    // The CloudToVTK functionObject from
    //
    //   lagrangian/intermediate/submodels/CloudFunctionObjects/
    //
    // stores file state on cloud OutputProperties itself,
    //
    // whereas the vtkCloud functionObject treats it like other
    // output and stores via the stateFunctionObject.
    // Since it uses VTP format, there is only a single file with all fields
    // - lookup by cloudName.

    const dictionary& cloudDict =
        geometryBase::parent_.mesh().lookupObject<IOdictionary>
        (
            cloudName_ + "OutputProperties"
        );

    if (cloudDict.found("cloudFunctionObject"))
    {
        const dictionary& foDict = cloudDict.subDict("cloudFunctionObject");
        if (foDict.found(functionObjectName_))
        {
            foDict.subDict(functionObjectName_)
                .readIfPresent("file", inputFileName_);
        }
    }
    else
    {
        inputFileName_ = getFileName("file", cloudName_);
    }


    stringOps::inplaceExpand(inputFileName_);
    if (inputFileName_.empty())
    {
        WarningInFunction
            << "Unable to find function object " << functionObjectName_
            << " output for field " << fieldName_
            << ". Cloud will not be processed"
            << endl;
        return;
    }

    vtkSmartPointer<vtkPolyData> dataset;

    if (inputFileName_.hasExt("vtp"))
    {
        auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
        reader->SetFileName(inputFileName_.c_str());
        reader->Update();

        dataset = reader->GetOutput();
    }
    else if (inputFileName_.hasExt("vtk"))
    {
        auto reader = vtkSmartPointer<vtkPolyDataReader>::New();
        reader->SetFileName(inputFileName_.c_str());
        reader->Update();

        dataset = reader->GetOutput();
    }
    else
    {
        // Invalid name - ignore
        inputFileName_.clear();
    }

    if (dataset)
    {
        auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();

        actor_->SetMapper(mapper);

        addGlyphs
        (
            position,
            fieldName_,
            colourFieldName_,
            maxGlyphLength_,
            dataset,
            actor_,
            renderer
        );

        renderer->AddActor(actor_);
    }
}


void Foam::functionObjects::runTimePostPro::functionObjectCloud::updateActors
(
    const scalar position
)
{
    actor_->GetProperty()->SetOpacity(opacity(position));

    vector pc = pointColour_->value(position);
    actor_->GetProperty()->SetColor(pc[0], pc[1], pc[2]);
}


bool Foam::functionObjects::runTimePostPro::functionObjectCloud::clear()
{
    if (functionObjectBase::clear())
    {
        if (inputFileName_.size() && Foam::rm(inputFileName_))
        {
            inputFileName_.clear();
            return true;
        }
    }

    return false;
}


// ************************************************************************* //
