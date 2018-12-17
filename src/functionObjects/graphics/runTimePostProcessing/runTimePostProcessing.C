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


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
static void addGeometryToScene
(
    PtrList<Type>& objects,
    const scalar position,
    vtkRenderer* renderer
)
{
    for (Type& obj : objects)
    {
        obj.addGeometryToScene(position, renderer);
    }
}


template<class Type>
static void updateActors(PtrList<Type>& objects, const scalar position)
{
    for (Type& obj : objects)
    {
        obj.updateActors(position);
    }
}


template<class Type>
static void cleanup(PtrList<Type>& objects)
{
    for (Type& obj : objects)
    {
        obj.clear();
    }
}

} // End namespace Foam


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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::runTimePostProcessing::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    Info<< type() << " " << name() << ": reading post-processing data" << endl;

    scene_.read(dict);

    const dictionary& outputDict = dict.subDict("output");
    outputDict.readEntry("name", output_.name_);
    outputDict.readEntry("width", output_.width_);
    outputDict.readEntry("height", output_.height_);

    readObjects(dict.subOrEmptyDict("points"), points_);
    readObjects(dict.subOrEmptyDict("lines"), lines_);
    readObjects(dict.subOrEmptyDict("surfaces"), surfaces_);

    const dictionary& textDict = dict.subDict("text");

    for (const entry& dEntry : textDict)
    {
        if (!dEntry.isDict())
        {
            FatalIOErrorInFunction(textDict)
                << textDict.dictName()
                << "text must be specified in dictionary format"
                << exit(FatalIOError);
        }

        const dictionary& objectDict = dEntry.dict();

        text_.append
        (
            new runTimePostPro::text
            (
                *this,
                objectDict,
                scene_.colours()
            )
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


    // Disable any floating point trapping
    // (some low-level rendering functionality does not like it)

    sigFpe::ignore sigFpeHandling; //<- disable in local scope

    // Initialise render window
    auto renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->OffScreenRenderingOn();
    renderWindow->SetSize(output_.width_, output_.height_);

    // Legacy rendering - was deprecated for 8.1.0
    #if (VTK_MAJOR_VERSION < 8) || \
        ((VTK_MAJOR_VERSION == 8) && (VTK_MINOR_VERSION < 2))
    renderWindow->SetAAFrames(10);
    #endif
    renderWindow->SetAlphaBitPlanes(true);
    renderWindow->SetMultiSamples(0);
//    renderWindow->PolygonSmoothingOn();

    auto renderer = vtkSmartPointer<vtkRenderer>::New();
    scene_.initialise(renderer, output_.name_);

    renderWindow->AddRenderer(renderer);


    addGeometryToScene(points_, 0, renderer);
    addGeometryToScene(lines_, 0, renderer);
    addGeometryToScene(surfaces_, 0, renderer);
    addGeometryToScene(text_, 0, renderer);

    while (scene_.loop(renderer))
    {
        const scalar position = scene_.position();

        updateActors(text_, position);
        updateActors(points_, position);
        updateActors(lines_, position);
        updateActors(surfaces_, position);
    }

    // Cleanup
    cleanup(text_);
    cleanup(points_);
    cleanup(lines_);
    cleanup(surfaces_);


    // Instead of relying on the destructor, manually restore the previous
    // SIGFPE state.
    // This is only to avoid compiler complaints about unused variables.

    sigFpeHandling.restore();

    return true;
}


// ************************************************************************* //
