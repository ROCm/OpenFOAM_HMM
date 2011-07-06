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

#include "waterFluxController.H"
#include "IFstream.H"
#include "graph.H"
#include "polyatomicCloud.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(waterFluxController, 0);

defineRunTimeSelectionTable(waterFluxController, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
waterFluxController::waterFluxController
(
    Time& t,
    polyatomicCloud& cloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(cloud.mesh())),
    cloud_(cloud),
    rndGen_(clock::getTime()),
    controllerDict_(dict.subDict("controllerProperties")),
    timeDict_(controllerDict_.subDict("timeProperties")),
    time_(t, timeDict_),
    regionName_(controllerDict_.lookup("zoneName")),
    regionId_(-1),
    zoneSurfaceArea_(0.0),
    internalFaces_(),
    processorFaces_(),
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
    const faceZoneMesh& faceZones = mesh_.faceZones();
    regionId_ = faceZones.findZoneID(regionName_);

    if (regionId_ == -1)
    {
        FatalErrorIn("waterFluxController::waterFluxController()")
        << "Cannot find region (faceZone): " << regionName_ << nl << "in: "
        << t.time().system()/"controllersDict"
        << exit(FatalError);
    }

    control_ = Switch(controllerDict_.lookup("controlSwitch"));
    readStateFromFile_ = Switch(controllerDict_.lookup("readStateFromFile"));

    setFacesInfo();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<waterFluxController> waterFluxController::New
(
    Time& t,
    polyatomicCloud& cloud,
    const dictionary& dict
)
{
    word waterFluxControllerName
    (
        dict.lookup("fluxControllerModel")
    );

    Info<< "Selecting fluxController "
         << waterFluxControllerName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(waterFluxControllerName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "waterFluxController::New(const dictionary&) : " << endl
            << "    unknown waterFluxController type "
            << waterFluxControllerName
            << ", constructor not in hash table" << endl << endl
            << "    Valid injector types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<waterFluxController>
        (
                cstrIter()(t, cloud, dict)
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

waterFluxController::~waterFluxController()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// void waterFluxController::updateTime()
// {
//     time_++;
//
//     const scalar& t = time_.time().timeOutputValue();
//
//     if ((t - initialTime_) < timePeriod_)
//     {
//         time_.controlTimeInterval().endTime() = false;
// //         control_ = false;
//     }
//     else
//     {
// //         control_ = true;
//     }
// }


void waterFluxController::setFacesInfo()
{
    const labelList& faces = controlZone();

    if (Pstream::parRun())
    {
        DynamicList<label> processorFaces(0);

        forAll(mesh_.boundaryMesh(), patchI)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[patchI];

            if (isA<processorPolyPatch>(patch))
            {
                for (label p = 0; p < patch.size(); p++)
                {
                    label patchFaceI = p + patch.start();
                    label faceId = findIndex (faces, patchFaceI);

                    if (faceId != -1)
                    {
                        processorFaces.append(patchFaceI);
                    }
                }
            }
        }

        processorFaces.shrink();

        processorFaces_.setSize(processorFaces.size(), -1);

        forAll(processorFaces, f)
        {
            processorFaces_[f] = processorFaces[f];
        }

        label nInternalFaces = faces.size() - processorFaces.size();
        internalFaces_.setSize(nInternalFaces, -1);

        label counter = 0;

        forAll(faces, f)
        {
            const label& faceI = faces[f];

            if (findIndex(processorFaces, faceI) == -1)
            {
                internalFaces_[counter] = faceI;
                counter++;
            }
        }

//         Pout << "processorFaces: " << processorFaces_ << endl;
//         Pout << "internalFaces: " << internalFaces_ << endl;

        forAll(internalFaces_, f)
        {
            const label& faceI = internalFaces_[f];
            zoneSurfaceArea_ += mag(mesh_.faceAreas()[faceI]);
        }

        // faces on a zone located on a processor cut belong to both processors
        // (hence the 0.5)

        forAll(processorFaces_, f)
        {
            const label& faceI = processorFaces_[f];
            zoneSurfaceArea_ += 0.5*mag(mesh_.faceAreas()[faceI]);
        }


        if (Pstream::parRun())
        {
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if (p != Pstream::myProcNo())
                {
                    const int proc = p;
                    {
                        OPstream toNeighbour(Pstream::blocking, proc);
                        toNeighbour << zoneSurfaceArea_;
                    }
                }
            }

            //- receiving
            for (int p = 0; p < Pstream::nProcs(); p++)
            {
                if (p != Pstream::myProcNo())
                {
                    scalar zoneSurfaceAreaProc;

                    const int proc = p;
                    {
                        IPstream fromNeighbour(Pstream::blocking, proc);
                        fromNeighbour >> zoneSurfaceAreaProc;
                    }

                    zoneSurfaceArea_ += zoneSurfaceAreaProc;
                }
            }
        }
    }
    else
    {
        forAll(faces, f)
        {
            const label& faceI = faces[f];

            zoneSurfaceArea_ += mag(mesh_.faceAreas()[faceI]);
        }
    }
}



void waterFluxController::updateTime()
{
    time_++;

//     const scalar& t = time_.time().timeOutputValue();
//
//     if ((t - initialTime_) < timePeriod_)
//     {
//         time_.controlTimeInterval().endTime() = false;
// //         control_ = false;
//     }
//     else
//     {
// //         control_ = true;
//     }
}

void waterFluxController::updateFluxControllerProperties
(
    const dictionary& newDict
)
{
    controllerDict_ = newDict.subDict("controllerProperties");

    //- you can reset the controlling zone from here. This essentially
    //  means that the coupling zone can infact move arbitrarily. To make
    //  this happen we probably need to devise a technique for automatically
    //  changing the cellZone else where, and then calling this function to
    //  reset the controlling zone in which the controller operates in.

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
}

const labelList& waterFluxController::controlZone() const
{
    return mesh_.faceZones()[regionId_];
}

label waterFluxController::isFaceOnControlZone(const label& faceI)
{
    const label f = findIndex(controlZone(), faceI);

    return f;
}

const word& waterFluxController::regionName() const
{
    return regionName_;
}

const scalar& waterFluxController::density() const
{
    return density_;
}

scalar& waterFluxController::density()
{
    return density_;
}

const vector& waterFluxController::velocity() const
{
    return velocity_;
}

vector& waterFluxController::velocity()
{
    return velocity_;
}

const scalar& waterFluxController::temperature() const
{
    return temperature_;
}

scalar& waterFluxController::temperature()
{
    return temperature_;
}

const scalar& waterFluxController::pressure() const
{
    return pressure_;
}

scalar& waterFluxController::pressure()
{
    return pressure_;
}

const tensor& waterFluxController::strainRate() const
{
    return strainRate_;
}

tensor& waterFluxController::strainRate()
{
    return strainRate_;
}

const vector& waterFluxController::tempGradient() const
{
    return tempGradient_;
}

vector& waterFluxController::tempGradient()
{
    return tempGradient_;
}


const scalarField& waterFluxController::densityField() const
{
    return densities_;
}

scalarField& waterFluxController::densityField()
{
    return densities_;
}

const vectorField& waterFluxController::velocityField() const
{
    return velocities_;
}
vectorField& waterFluxController::velocityField()
{
    return velocities_;
}

const scalarField& waterFluxController::temperatureField() const
{
    return temperatures_;
}

scalarField& waterFluxController::temperatureField()
{
    return temperatures_;
}

const scalarField& waterFluxController::pressureField() const
{
    return pressures_;
}

scalarField& waterFluxController::pressureField()
{
    return pressures_;
}


const bool& waterFluxController::singleValueController() const
{
    return singleValueController_;
}

bool& waterFluxController::singleValueController()
{
    return singleValueController_;
}

const bool& waterFluxController::fieldController() const
{
    return fieldController_;
}

bool& waterFluxController::fieldController()
{
    return fieldController_;
}


const bool& waterFluxController::writeInTimeDir() const
{
    return writeInTimeDir_;
}

const bool& waterFluxController::writeInCase() const
{
    return writeInCase_;
}


// const scalar waterFluxController::avReqDensity() const
// {
//     scalar totalDensity = 0.0;
//
//     forAll(densities_, c)
//     {
//         totalDensity += densities_[c];
//     }
//
//     if (cells_.size() > 0)
//     {
//         totalDensity /= scalar(cells_.size());
//     }
//
//     return totalDensity;
// }
//
// const vector waterFluxController::avReqVelocity() const
// {
//     vector totalVel = vector::zero;
//
//     forAll(velocities_, c)
//     {
//         totalVel += velocities_[c];
//     }
//
//     if (cells_.size() > 0)
//     {
//         totalVel /= scalar(cells_.size());
//     }
//
//     return totalVel;
// }
//
// const scalar waterFluxController::avReqTemperature() const
// {
//     scalar totalTemp = 0.0;
//
//     forAll(densities_, c)
//     {
//         totalTemp += temperatures_[c];
//     }
//
//     if (cells_.size() > 0)
//     {
//         totalTemp /= scalar(cells_.size());
//     }
//
//     return totalTemp;
// }



} // End namespace Foam

// ************************************************************************* //
