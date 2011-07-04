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

#include "stateController.H"
#include "IFstream.H"
#include "polyatomicCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(stateController, 0);

defineRunTimeSelectionTable(stateController, dictionary);

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::stateController::stateController
(
    polyatomicCloud& cloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(cloud.mesh())),
    cloud_(cloud),
    rndGen_(clock::getTime()),
    controllerDict_(dict.subDict("controllerProperties")),
    timeDict_(controllerDict_.subDict("timeProperties")),
    time_(mesh_.time(), timeDict_),
    timePeriod_(readScalar(timeDict_.lookup("initialTimePeriod"))), //temp
    initialTime_(time_.time().startTime().value()),
    regionName_(controllerDict_.lookup("zoneName")),
    regionId_(-1),
    control_(true),
    readStateFromFile_(true),
    singleValueController_(false),
    density_(0.0),
    velocity_(vector::zero),
    temperature_(0.0),
    pressure_(0.0),
    strainRate_(tensor::zero),
    tempGradient_(vector::zero),
    fieldController_(false),
    densities_(),
    velocities_(),
    temperatures_(),
    pressures_(),
    writeInTimeDir_(true),
    writeInCase_(true)
{
    const cellZoneMesh& cellZones = mesh_.cellZones();

    regionId_ = cellZones.findZoneID(regionName_);

    if (regionId_ == -1)
    {
        FatalErrorIn("stateController::stateController()")
            << "Cannot find region: " << regionName_ << nl << "in: "
            << time_.time().system()/"controllersDict"
            << exit(FatalError);
    }

    control_ = Switch(controllerDict_.lookup("controlSwitch"));

    readStateFromFile_ = Switch(controllerDict_.lookup("readStateFromFile"));

    const scalar& avTimeInterval = time_.averageTimeInterval().deltaT();

    if ((timePeriod_ < avTimeInterval) && (timePeriod_ > 0.0))
    {
        timePeriod_ = avTimeInterval;
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::stateController> Foam::stateController::New
(
    polyatomicCloud& cloud,
    const dictionary& dict
)
{
    word stateControllerName
    (
        dict.lookup("stateControllerModel")
    );

    Info<< "Selecting stateController "
         << stateControllerName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(stateControllerName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "stateController::New(const dictionary&) : " << endl
            << "    unknown stateController type "
            << stateControllerName
            << ", constructor not in hash table" << endl << endl
            << "    Valid types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<stateController>
    (
        cstrIter()(cloud, dict)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::stateController::~stateController()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::stateController::updateTime()
{
    time_++;

    const scalar& t = time_.time().timeOutputValue();

    if ((t - initialTime_) < timePeriod_)
    {
        time_.controlTimeInterval().endTime() = false;
    }
}


void Foam::stateController::updateStateControllerProperties
(
    const dictionary& newDict
)
{
    controllerDict_ = newDict.subDict("controllerProperties");

    if (controllerDict_.found("controlSwitch"))
    {
        control_ = Switch(controllerDict_.lookup("controlSwitch"));
    }

    if (controllerDict_.found("readStateFromFile"))
    {
        readStateFromFile_ = Switch
        (
            controllerDict_.lookup("readStateFromFile")
        );
    }

    timeDict_ = controllerDict_.subDict("timeProperties");

    if (timeDict_.found("resetAtOutput"))
    {
        time_.resetFieldsAtOutput() = Switch(timeDict_.lookup("resetAtOutput"));
    }
}


const Foam::labelList& Foam::stateController::controlZone() const
{
    return mesh_.cellZones()[regionId_];
}

const Foam::word& Foam::stateController::regionName() const
{
    return regionName_;
}


Foam::scalar Foam::stateController::density() const
{
    return density_;
}


Foam::scalar& Foam::stateController::density()
{
    return density_;
}


const Foam::vector& Foam::stateController::velocity() const
{
    return velocity_;
}


Foam::vector& Foam::stateController::velocity()
{
    return velocity_;
}


Foam::scalar Foam::stateController::temperature() const
{
    return temperature_;
}


Foam::scalar& Foam::stateController::temperature()
{
    return temperature_;
}

const Foam::scalar& Foam::stateController::pressure() const
{
    return pressure_;
}


Foam::scalar& Foam::stateController::pressure()
{
    return pressure_;
}


const Foam::tensor& Foam::stateController::strainRate() const
{
    return strainRate_;
}


Foam::tensor& Foam::stateController::strainRate()
{
    return strainRate_;
}


const Foam::vector& Foam::stateController::tempGradient() const
{
    return tempGradient_;
}


Foam::vector& Foam::stateController::tempGradient()
{
    return tempGradient_;
}


const Foam::scalarField& Foam::stateController::densityField() const
{
    return densities_;
}

Foam::scalarField& Foam::stateController::densityField()
{
    return densities_;
}


const Foam::vectorField& Foam::stateController::velocityField() const
{
    return velocities_;
}


Foam::vectorField& Foam::stateController::velocityField()
{
    return velocities_;
}


const Foam::scalarField& Foam::stateController::temperatureField() const
{
    return temperatures_;
}


Foam::scalarField& Foam::stateController::temperatureField()
{
    return temperatures_;
}


const Foam::scalarField& Foam::stateController::pressureField() const
{
    return pressures_;
}


Foam::scalarField& Foam::stateController::pressureField()
{
    return pressures_;
}


bool Foam::stateController::singleValueController() const
{
    return singleValueController_;
}


bool& Foam::stateController::singleValueController()
{
    return singleValueController_;
}


bool Foam::stateController::fieldController() const
{
    return fieldController_;
}


bool& Foam::stateController::fieldController()
{
    return fieldController_;
}


bool Foam::stateController::writeInTimeDir() const
{
    return writeInTimeDir_;
}


bool Foam::stateController::writeInCase() const
{
    return writeInCase_;
}


Foam::scalar Foam::stateController::avReqDensity()
{
    scalar totalDensity = 0.0;

    if (singleValueController_)
    {
        totalDensity = density_;
    }
    else if (fieldController_)
    {
        label controlCells = controlZone().size();

        forAll(densities_, c)
        {
            totalDensity += densities_[c];
        }

        if (Pstream::parRun())
        {
            reduce(totalDensity, sumOp<scalar>());

            reduce(controlCells, sumOp<label>());
        }

        if (controlCells > 0)
        {
            totalDensity /= scalar(controlCells);
        }
    }

    return totalDensity;
}


Foam::vector Foam::stateController::avReqVelocity()
{
    vector totalVel = vector::zero;

    if (singleValueController_)
    {
        totalVel = velocity_;
    }
    else if (fieldController_)
    {
        label controlCells = controlZone().size();

        forAll(velocities_, c)
        {
            totalVel += velocities_[c];
        }

        if (Pstream::parRun())
        {
            reduce(totalVel, sumOp<vector>());

            reduce(controlCells, sumOp<label>());
        }

        if (controlCells > 0)
        {
            totalVel /= scalar(controlCells);
        }
    }

    return totalVel;
}


Foam::scalar Foam::stateController::avReqTemperature()
{
    scalar totalTemp = 0.0;

    if (singleValueController_)
    {
        totalTemp = temperature_;
    }
    else if (fieldController_)
    {
        label controlCells = controlZone().size();

        forAll(temperatures_, c)
        {
            totalTemp += temperatures_[c];
        }

        if (Pstream::parRun())
        {
            reduce(totalTemp, sumOp<scalar>());

            reduce(controlCells, sumOp<label>());
        }

        if (controlCells > 0)
        {
            totalTemp /= scalar(controlCells);
        }
    }

    return totalTemp;
}


Foam::scalar Foam::stateController::avReqPressure()
{
    scalar totalPressure = 0.0;

    if (singleValueController_)
    {
        totalPressure = pressure_;
    }
    else if (fieldController_)
    {
        label controlCells = controlZone().size();

        forAll(pressures_, c)
        {
            totalPressure += pressures_[c];
        }

        if (Pstream::parRun())
        {
            reduce(totalPressure, sumOp<scalar>());

            reduce(controlCells, sumOp<label>());
        }

        if (controlCells > 0)
        {
            totalPressure /= scalar(controlCells);
        }
    }

    return totalPressure;
}


// ************************************************************************* //
