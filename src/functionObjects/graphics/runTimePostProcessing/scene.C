/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2017 OpenCFD Ltd.
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
#include "scene.H"
#include "Constant.H"

// VTK includes
#include "vtkCamera.h"
#include "vtkCubeSource.h"
#include "vtkLightKit.h"
#include "vtkPolyDataMapper.h"
#include "vtkPNGWriter.h"
#include "vtkRenderer.h"
#include "vtkRendererCollection.h"
#include "vtkRenderWindow.h"
#include "vtkWindowToImageFilter.h"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::runTimePostPro::scene::readCamera
(
    const dictionary& dict
)
{
    if (dict.readIfPresent("nFrameTotal", nFrameTotal_))
    {
        if (nFrameTotal_ < 1)
        {
            FatalIOErrorInFunction(dict)
                << "nFrameTotal must be 1 or greater"
                << exit(FatalIOError);
        }
    }

    if (dict.readIfPresent("startPosition", startPosition_))
    {
        if ((startPosition_ < 0) || (startPosition_ > 1))
        {
            FatalIOErrorInFunction(dict)
                << "startPosition must be in the range 0-1"
                << exit(FatalIOError);
        }
        else
        {
            position_ = startPosition_;
        }
    }

    if (nFrameTotal_ > 1)
    {
        scalar endPosition = dict.lookupOrDefault<scalar>("endPosition", 1);
        if ((endPosition < 0) || (endPosition > 1))
        {
            FatalIOErrorInFunction(dict)
                << "endPosition must be in the range 0-1"
                << exit(FatalIOError);
        }
        dPosition_ = (endPosition - startPosition_)/scalar(nFrameTotal_ - 1);
    }

    cameraPosition_ = Function1<vector>::New("position", dict);
    cameraFocalPoint_ = Function1<point>::New("focalPoint", dict);
    cameraUp_ = Function1<vector>::New("up", dict);

    dict.readIfPresent("clipBox", clipBox_);
    dict.lookup("parallelProjection") >> parallelProjection_;
    if (!parallelProjection_)
    {
        if (dict.found("viewAngle"))
        {
            cameraViewAngle_ = Function1<scalar>::New("viewAngle", dict);
        }
        else
        {
            cameraViewAngle_.reset
            (
                new Function1Types::Constant<scalar>("viewAngle", 35.0)
            );
        }
    }

    if (dict.found("zoom"))
    {
        cameraZoom_ = Function1<scalar>::New("zoom", dict);
    }
    else
    {
        cameraZoom_.reset
        (
            new Function1Types::Constant<scalar>("zoom", 1.0)
        );
    }
}


void Foam::functionObjects::runTimePostPro::scene::readColours
(
    const dictionary& dict
)
{
    const wordList colours = dict.toc();
    forAll(colours, i)
    {
        const word& c = colours[i];
        colours_.insert(c, Function1<vector>::New(c, dict).ptr());
    }
}


void Foam::functionObjects::runTimePostPro::scene::setActorVisibility
(
    vtkRenderer* renderer,
    const bool visible
) const
{
    vtkActorCollection *actors = renderer->GetActors();
    for (int i = 0; i < actors->GetNumberOfItems(); ++i)
    {
        vtkActor *actor = vtkActor::SafeDownCast(actors->GetItemAsObject(i));
        actor->SetVisibility(visible);
    }
}


void Foam::functionObjects::runTimePostPro::scene::initialise
(
    vtkRenderer* renderer,
    const word& outputName
)
{
    currentFrameI_ = 0;
    position_ = startPosition_;

    outputName_ = outputName;

    // Set the background
    const vector backgroundColour = colours_["background"]->value(position_);
    renderer->SetBackground
    (
        backgroundColour.x(),
        backgroundColour.y(),
        backgroundColour.z()
    );

    // Apply gradient background if "background2" defined
    if (colours_.found("background2"))
    {
        renderer->GradientBackgroundOn();
        vector backgroundColour2 = colours_["background2"]->value(position_);

        renderer->SetBackground2
        (
            backgroundColour2.x(),
            backgroundColour2.y(),
            backgroundColour2.z()
        );
    }

    // Depth peeling
    renderer->SetUseDepthPeeling(true);
    renderer->SetMaximumNumberOfPeels(4);
    renderer->SetOcclusionRatio(0);

    // Set the camera
    auto camera = vtkSmartPointer<vtkCamera>::New();
    camera->SetParallelProjection(parallelProjection_);
    renderer->SetActiveCamera(camera);

    // Add the lights
    auto lightKit = vtkSmartPointer<vtkLightKit>::New();
    lightKit->AddLightsToRenderer(renderer);

    if (!clipBox_.empty())
    {
        const point& min = clipBox_.min();
        const point& max = clipBox_.max();
        auto clipBox = vtkSmartPointer<vtkCubeSource>::New();
        clipBox->SetXLength(max.x() - min.x());
        clipBox->SetYLength(max.y() - min.y());
        clipBox->SetZLength(max.z() - min.z());
        clipBox->SetCenter
        (
            min.x() + 0.5*(max.x() - min.x()),
            min.y() + 0.5*(max.y() - min.y()),
            min.z() + 0.5*(max.z() - min.z())
        );
        auto clipMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        clipMapper->SetInputConnection(clipBox->GetOutputPort());

        clipBoxActor_ = vtkSmartPointer<vtkActor>::New();
        clipBoxActor_->SetMapper(clipMapper);
        clipBoxActor_->VisibilityOff();
        renderer->AddActor(clipBoxActor_);
    }
}


void Foam::functionObjects::runTimePostPro::scene::setCamera
(
    vtkRenderer* renderer
) const
{
    vtkCamera* camera = renderer->GetActiveCamera();

    if (parallelProjection_)
    {
        // Restore parallel scale to allow application of zoom (later)
        camera->SetParallelScale(1);
    }
    else
    {
        // Restore viewAngle (it might be reset by clipping)
        camera->SetViewAngle(cameraViewAngle_->value(position_));
    }

    const vector up = cameraUp_->value(position_);
    const vector pos = cameraPosition_->value(position_);
    const point focalPoint = cameraFocalPoint_->value(position_);
    const scalar zoom = cameraZoom_->value(position_);

    camera->SetViewUp(up.x(), up.y(), up.z());
    camera->SetPosition(pos.x(), pos.y(), pos.z());
    camera->SetFocalPoint(focalPoint.x(), focalPoint.y(), focalPoint.z());


    // Apply clipping if required
    // Note: possible optimisation - if the camera is static, this only needs
    //       to be done once on initialisation
    if (!clipBox_.empty())
    {
        setActorVisibility(renderer, false);
        clipBoxActor_->VisibilityOn();

        // Call ResetCamera() to fit clip box in view
        renderer->ResetCamera();

        setActorVisibility(renderer, true);
        clipBoxActor_->VisibilityOff();
    }

    // Zoom applied after all other operations
    camera->Zoom(zoom);

    camera->Modified();
}


Foam::string
Foam::functionObjects::runTimePostPro::scene::frameIndexStr() const
{
    string str = Foam::name(currentFrameI_);
    str.insert(0, 4 - str.length(), '0');

    return str;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::scene::scene
(
    const objectRegistry& obr,
    const word& name
)
:
    obr_(obr),
    name_(name),
    colours_(),
    cameraPosition_(nullptr),
    cameraFocalPoint_(nullptr),
    cameraUp_(nullptr),
    cameraViewAngle_(nullptr),
    cameraZoom_(nullptr),
    clipBox_(boundBox::invertedBox),
    clipBoxActor_(),
    parallelProjection_(true),
    nFrameTotal_(1),
    startPosition_(0),
    position_(0),
    dPosition_(0),
    currentFrameI_(0),
    outputName_("unknown")
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::scene::~scene()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::HashPtrTable<Foam::Function1<Foam::vector>, Foam::word>&
Foam::functionObjects::runTimePostPro::scene::colours() const
{
    return colours_;
}


Foam::label Foam::functionObjects::runTimePostPro::scene::frameIndex() const
{
    return currentFrameI_;
}


Foam::scalar Foam::functionObjects::runTimePostPro::scene::position() const
{
    return position_;
}


void Foam::functionObjects::runTimePostPro::scene::read
(
    const dictionary& dict
)
{
    readCamera(dict.subDict("camera"));
    readColours(dict.subDict("colours"));
}


bool Foam::functionObjects::runTimePostPro::scene::loop(vtkRenderer* renderer)
{
    static bool initialised = false;
    setCamera(renderer);

    if (!initialised)
    {
        initialised = true;
        return true;
    }

    // Ensure that all objects can be seen without clipping
    // Note: can only be done after all objects have been added!
    renderer->ResetCameraClippingRange();

    // Save image from last iteration
    saveImage(renderer->GetRenderWindow());

    currentFrameI_++;

    position_ = startPosition_ + currentFrameI_*dPosition_;

    if (currentFrameI_ < nFrameTotal_)
    {
        return true;
    }
    else
    {
        initialised = false;
        return false;
    }
}


void Foam::functionObjects::runTimePostPro::scene::saveImage
(
    vtkRenderWindow* renderWindow
) const
{
    if (!renderWindow)
    {
        return;
    }

    const Time& runTime = obr_.time();

    fileName prefix(Pstream::parRun() ?
        runTime.path()/".."/"postProcessing"/name_/obr_.time().timeName() :
        runTime.path()/"postProcessing"/name_/obr_.time().timeName());

    mkDir(prefix);

    renderWindow->Render();

    // Set up off-screen rendering
    auto windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();

    windowToImageFilter->SetInput(renderWindow);

    //// Add alpha channel for transparency
    // windowToImageFilter->SetInputBufferTypeToRGBA();
    windowToImageFilter->SetInputBufferTypeToRGB();

//    windowToImageFilter->ReadFrontBufferOff();
    windowToImageFilter->Update();

    // Save the image
    auto writer = vtkSmartPointer<vtkPNGWriter>::New();
    fileName fName(prefix/outputName_ + '.' + frameIndexStr() + ".png");
    writer->SetFileName(fName.c_str());
    writer->SetInputConnection(windowToImageFilter->GetOutputPort());

    Info<< "    Generating image: " << fName << endl;

    writer->Write();
}


// ************************************************************************* //
