/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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

#include "controllers.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::controllers::controllers
(
    const polyMesh& mesh
)
:
    time_(mesh.time()),
    controllersDict_
    (
        IOobject
        (
            "controllersDict",
            time_.system(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    stateControllersList_(),
    sCNames_(),
    sCIds_(),
    sCFixedPathNames_(),
    stateControllers_(),
    fluxControllersList_(),
    fCNames_(),
    fCIds_(),
    fCFixedPathNames_(),
    fluxControllers_()
{}


Foam::controllers::controllers
(
    const polyMesh& mesh,
    polyatomicCloud& cloud
)
:
    time_(mesh.time()),
    controllersDict_
    (
        IOobject
        (
            "controllersDict",
            time_.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    stateControllersList_(controllersDict_.lookup("stateControllers")),
    sCNames_(stateControllersList_.size()),
    sCIds_(stateControllersList_.size()),
    sCFixedPathNames_(stateControllersList_.size()),
    stateControllers_(stateControllersList_.size()),
    fluxControllersList_(controllersDict_.lookup("fluxControllers")),
    fCNames_(fluxControllersList_.size()),
    fCIds_(fluxControllersList_.size()),
    fCFixedPathNames_(fluxControllersList_.size()),
    fluxControllers_(fluxControllersList_.size())
{

    Info << nl << "Creating controllers" << nl << endl;

    // state controllers

    if (!stateControllers_.empty())
    {
        forAll(stateControllers_, sC)
        {
            const entry& controllersI = stateControllersList_[sC];

            const dictionary& controllersIDict = controllersI.dict();

            stateControllers_[sC] = autoPtr<stateController>
            (
                stateController::New(time_, cloud, controllersIDict)
            );

            sCNames_[sC] = stateControllers_[sC]->type();

            sCIds_[sC] = sC;
        }
    }

    //- flux controllers

    if (!fluxControllers_.empty())
    {
        forAll(fluxControllers_, fC)
        {
            const entry& controllersI = fluxControllersList_[fC];

            const dictionary& controllersIDict = controllersI.dict();

            fluxControllers_[fC] = autoPtr<fluxController>
            (
                fluxController::New(time_, cloud, controllersIDict)
            );

            fCNames_[fC] = fluxControllers_[fC]->type();
            fCIds_[fC] = fC;
        }
    }

    // creating directories for state controllers
    if (!nStateControllers_.empty())
    {
        // case/controllers
        fileName controllersPath(time_.path()/"controllers");

        if (!isDir(controllersPath))
        {
            mkDir(controllersPath);
        }

        // case/controllers/<cloudName>
        fileName controllersPath(controllersPath/cloud.name());

        if (!isDir(controllersPath))
        {
            mkDir(controllersPath);
        }

        // case/controllers/<cloudName>/stateControllers
        fileName stateControllersPath(controllersPath/"stateControllers");

        if (!isDir(stateControllersPath))
        {
            mkDir(stateControllersPath);
        }

        forAll(stateControllers_, sC)
        {
            if (stateControllers_[sC]->writeInCase())
            {
                // case/controllers/<cloudName>/
                //     stateControllers/<stateControllerModel>
                fileName stateControllerPath(stateControllersPath/sCNames_[sC]);

                if (!isDir(stateControllerPath))
                {
                    mkDir(stateControllerPath);
                }

                const word& regionName = stateControllers_[sC]->regionName();

                // case/controllers/<cloudName>/
                //     stateControllers/<stateControllerModel>/<cellZoneName>
                fileName zonePath(stateControllerPath/regionName);

                if (!isDir(zonePath))
                {
                    mkDir(zonePath);
                }

                sCFixedPathNames_[sC] = zonePath;
            }
        }
    }

    // creating directories for flux controllers
    if (nFluxControllers_ > 0)
    {
        // case/controllers
        fileName controllersPath(time_.path()/"controllers");

        if ( !isDir(controllersPath) )
        {
            mkDir(controllersPath);
        }

        // case/controllers/<cloudName>
        fileName controllersPath(time_.path()/cloud.name());

        if ( !isDir(controllersPath) )
        {
            mkDir(controllersPath);
        }

        // case/controllers/<cloudName>/fluxControllers
        fileName fluxControllersPath(controllersPath/"fluxControllers");

        if (!isDir(fluxControllersPath))
        {
            mkDir(fluxControllersPath);
        }

        forAll(fluxControllers_, fC)
        {
            if (fluxControllers_[fC]->writeInCase())
            {
                // case/controllers/<cloudName>/
                //    fluxControllers/<fluxControllerModel>
                fileName fluxControllerPath(fluxControllersPath/fCNames_[fC]);

                if (!isDir(fluxControllerPath))
                {
                    mkDir(fluxControllerPath);
                }

                const word& regionName = fluxControllers_[fC]->regionName();

                // case/controllers/<cloudName>/
                //    fluxControllers/<fluxControllerModel>/<faceZoneName>
                fileName zonePath(fluxControllerPath/regionName);

                if (!isDir(zonePath))
                {
                    mkDir(zonePath);
                }

                fCFixedPathNames_[fC] = zonePath;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

controllers::~controllers()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void controllers::initialConfig()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->initialConfiguration();
    }

    forAll(fluxControllers_, fC)
    {
        fluxControllers_[fC]->initialConfiguration();
    }
}


void controllers::updateTimeInfo()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->updateTime();
    }

    forAll(fluxControllers_, fC)
    {
        fluxControllers_[fC]->updateTime();
    }
}


void controllers::controlState()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->controlMols();
    }
}


void controllers::controlVelocitiesI()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->controlMolsBeg();
    }
}

void controllers::controlVelocitiesII()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->controlMolsEnd();
    }
}


void  controllers::controlPriorToForces()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->controlBeforeForces();
    }
}


void controllers::calculateStateProps()
{
    forAll(stateControllers_, sC)
    {
        stateControllers_[sC]->calculateProperties();
    }

    forAll(fluxControllers_, fC)
    {
        fluxControllers_[fC]->calculateProperties();
    }
}


void controllers::outputStateResults()
{
    const Time& runTime = time_;

    if (runTime.outputTime())
    {
        // creating a set of directories in the current time directory
        {
            List<fileName> timePathNames(sCFixedPathNames_.size());

            if (nStateControllers_ > 0)
            {
                if (Pstream::master())
                {
                    // case/<timeDir>/uniform
                    fileName uniformTimePath
                    (
                        runTime.path()/runTime.timeName()/"uniform"
                    );

                    if (!isDir(uniformTimePath))
                    {
                        mkDir(uniformTimePath);
                    }

                    if (!stateControllers_.empty())
                    {
                        // case/<timeDir>/uniform/controllers
                        fileName controllersTimePath
                        (
                            uniformTimePath/"controllers"
                        );

                        if (!isDir(controllersTimePath))
                        {
                            mkDir(controllersTimePath);
                        }

                        // case/<timeDir>/uniform/controllers/<cloudName>
                        fileName cloudTimePath
                        (
                            controllersTimePath/cloud.name()
                        );

                        if (!isDir(cloudTimePath))
                        {
                            mkDir(cloudTimePath);
                        }

                        // case/<timeDir>/uniform/controllers/<cloudName>/
                        fileName stateControllersTimePath
                        (
                            cloudTimePath/"stateControllers"
                        );

                        if (!isDir(stateControllersTimePath))
                        {
                            mkDir(stateControllersTimePath);
                        }

                        forAll(stateControllers_, sC)
                        {
                            if (stateControllers_[sC]->writeInTimeDir())
                            {
                                // case/<timeDir>/uniform/controllers/
                                //     <cloudName>/<stateControllerModel>
                                fileName sCTimePath
                                (
                                    stateControllersTimePath/sCNames_[sC]
                                );

                                if (!isDir(sCTimePath))
                                {
                                    mkDir(sCTimePath);
                                }

                                // Creating directory for different zones but
                                // of the same model
                                const word& regionName =
                                    stateControllers_[sC]->regionName();

                                // case/<timeDir>/uniform/controllers/
                                //    <cloudName>/<stateControllerModel>/
                                //    <cellZoneName>
                                fileName zoneTimePath(sCTimePath/regionName);

                                if (!isDir(zoneTimePath))
                                {
                                    mkDir(zoneTimePath);
                                }

                                timePathNames[sC] = zoneTimePath;
                            }
                        }
                    }
                }
            }

            // write out data
            forAll(stateControllers_, sC)
            {
                stateControllers_[sC]->output
                (
                    sCFixedPathNames_[sC],
                    timePathNames[sC]
                );
            }
        }

        {
            List<fileName> timePathNames(fCFixedPathNames_.size());

            if (nFluxControllers_ > 0)
            {
                if (Pstream::master())
                {
                    // case/<timeDir>/uniform
                    fileName uniformTimePath
                    (
                        runTime.path()/runTime.timeName()/"uniform"
                    );

                    if (!isDir(uniformTimePath))
                    {
                        mkDir(uniformTimePath);
                    }

                    if (!fluxControllers_.empty())
                    {
                        // case/<timeDir>/uniform/controllers
                        fileName controllersTimePath
                        (
                            uniformTimePath/"controllers"
                        );

                        if (!isDir(controllersTimePath))
                        {
                            mkDir(controllersTimePath);
                        }

                        // case/<timeDir>/uniform/controllers/<cloudName>
                        fileName cloudTimePath
                        (
                            controllersTimePath/cloud.name()
                        );

                        if (!isDir(cloudTimePath))
                        {
                            mkDir(cloudTimePath);
                        }

                        // case/<timeDir>/uniform/fluxControllers
                        fileName controllersTimePath
                        (
                            cloudTimePath/"fluxControllers"
                        );

                        if (!isDir(controllersTimePath))
                        {
                            mkDir(controllersTimePath);
                        }

                        forAll(fluxControllers_, fC)
                        {
                            if (stateControllers_[fC]->writeInTimeDir())
                            {
                                // case/<timeDir>/uniform/controllers/
                                //     <cloudName>/<fluxControllerModel>
                                fileName fCTimePath
                                (
                                    controllersTimePath/fCNames_[fC]
                                );

                                if (!isDir(fCTimePath))
                                {
                                    mkDir(fCTimePath);
                                }

                                const word& regionName =
                                    fluxControllers_[fC]->regionName();

                                // case/<timeDir>/uniform/controllers/
                                //     <cloudName>/<fluxControllerModel>/
                                //     <faceZoneName>
                                fileName zoneTimePath(fCTimePath/regionName);

                                if (!isDir(zoneTimePath))
                                {
                                    mkDir(zoneTimePath);
                                }

                                timePathNames[fC] = zoneTimePath;
                            }
                        }
                    }
                }
            }

            // write out data
            forAll(fluxControllers_, fC)
            {
                fluxControllers_[fC]->output
                (
                    fCFixedPathNames_[fC],
                    timePathNames[fC]
                );
            }
        }

        // Re-read dictionaries for modified properties (run-time selection)
        {
            stateControllersList_.clear();

            stateControllersList_ = controllersDict_.lookup("stateControllers");

            forAll(stateControllers_, sC)
            {
                const entry& controllersI = stateControllersList_[sC];
                const dictionary& controllersIDict = controllersI.dict();

                stateControllers_[sC]->updateProperties(controllersIDict);
            }
        }

        {
            fluxControllersList_.clear();

            fluxControllersList_ = controllersDict_.lookup("fluxControllers");

            forAll(fluxControllers_, fC)
            {
                const entry& controllersI = fluxControllersList_[fC];
                const dictionary& controllersIDict = controllersI.dict();

                fluxControllers_[fC]->updateProperties(controllersIDict);
            }
        }
    }
}


// ************************************************************************* //
