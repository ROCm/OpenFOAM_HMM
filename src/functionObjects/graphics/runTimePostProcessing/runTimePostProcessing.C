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
#include "runTimePostProcessing.H"
#include "dictionary.H"
#include "pointData.H"
#include "pathline.H"
#include "surface.H"
#include "text.H"
#include "Time.H"
#include "sigFpe.H"
#include "addToRunTimeSelectionTable.H"

// VTK includes
#include "vtkPolyDataMapper.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkSmartPointer.h"

#include "vtkLight.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(runTimePostProcessing, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        runTimePostProcessing,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostProcessing::runTimePostProcessing
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    scene_(runTime, name),
    points_(),
    lines_(),
    surfaces_(),
    text_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostProcessing::~runTimePostProcessing()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::runTimePostProcessing::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    Info<< type() << " " << name() << ": reading post-processing data" << endl;

    scene_.read(dict);

    const dictionary& outputDict = dict.subDict("output");
    outputDict.lookup("name") >> output_.name_;
    outputDict.lookup("width") >> output_.width_;
    outputDict.lookup("height") >> output_.height_;


    readObjects(dict.subOrEmptyDict("points"), points_);
    readObjects(dict.subOrEmptyDict("lines"), lines_);
    readObjects(dict.subOrEmptyDict("surfaces"), surfaces_);


    const dictionary& textDict = dict.subDict("text");
    forAllConstIter(dictionary, textDict, iter)
    {
        if (!iter().isDict())
        {
            FatalIOErrorInFunction(textDict)
                << "text must be specified in dictionary format"
                << exit(FatalIOError);
        }

        text_.append(new runTimePostPro::text
        (
            *this,
            iter().dict(),
            scene_.colours())
        );
    }

    return true;
}


bool Foam::functionObjects::runTimePostProcessing::execute()
{
    return true;
}


bool Foam::functionObjects::runTimePostProcessing::write()
{
    if (!Pstream::master())
    {
        return true;
    }

    Info<< type() << " " << name() <<  " output:" << nl
        << "    Constructing scene" << endl;

    // Unset any floating point trapping (some low-level rendering functionality
    // does not like it)
    sigFpe::unset(false);

    // Initialise render window
    auto renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->OffScreenRenderingOn();
    renderWindow->SetSize(output_.width_, output_.height_);
    renderWindow->SetAAFrames(10);
    renderWindow->SetAlphaBitPlanes(true);
    renderWindow->SetMultiSamples(0);
//    renderWindow->PolygonSmoothingOn();

    auto renderer = vtkSmartPointer<vtkRenderer>::New();
    scene_.initialise(renderer, output_.name_);

    renderWindow->AddRenderer(renderer);

    // Add the points
    forAll(points_, i)
    {
        points_[i].addGeometryToScene(0, renderer);
    }

    // Add the lines
    forAll(lines_, i)
    {
        lines_[i].addGeometryToScene(0, renderer);
    }

    // Add the surfaces
    forAll(surfaces_, i)
    {
        surfaces_[i].addGeometryToScene(0, renderer);
    }

    // Add the text
    forAll(text_, i)
    {
        text_[i].addGeometryToScene(0, renderer);
    }

    while (scene_.loop(renderer))
    {
        scalar position = scene_.position();

        // Update the text
        forAll(text_, i)
        {
            text_[i].updateActors(position);
        }

        // Update the points
        forAll(points_, i)
        {
            points_[i].updateActors(position);
        }

        // Update the lines
        forAll(lines_, i)
        {
            lines_[i].updateActors(position);
        }

        // Update the surfaces
        forAll(surfaces_, i)
        {
            surfaces_[i].updateActors(position);
        }
    }

    // Clean up
    forAll(text_, i)
    {
        text_[i].clear();
    }
    forAll(points_, i)
    {
        points_[i].clear();
    }
    forAll(lines_, i)
    {
        lines_[i].clear();
    }
    forAll(surfaces_, i)
    {
        surfaces_[i].clear();
    }

    // Reset any floating point trapping
    sigFpe::set(false);

    return true;
}


// ************************************************************************* //
