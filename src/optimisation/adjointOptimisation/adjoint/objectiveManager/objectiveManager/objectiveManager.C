/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2021 PCOpt/NTUA
    Copyright (C) 2013-2021 FOSS GP
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "objectiveManager.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(objectiveManager, 0);
defineRunTimeSelectionTable(objectiveManager, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

objectiveManager::objectiveManager
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
:
    regIOobject
    (
        IOobject
        (
            "objectiveManager" + adjointSolverName,
            mesh.time().system(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true  //register object
        )
    ),
    mesh_(mesh),
    dict_(dict),
    adjointSolverName_(adjointSolverName),
    primalSolverName_(primalSolverName),
    objectives_(0),
    weigthedObjectiveFile_(nullptr)
{
    // Construct objectives
    //~~~~~~~~~~~~~~~~~~~~~
    Info << "Constructing objective functions " << nl << endl;
    const word objectiveType = dict.get<word>("type");
    const dictionary& objectiveNamesDict(dict.subDict("objectiveNames"));
    wordList objectiveNames(objectiveNamesDict.toc());
    objectives_.setSize(objectiveNames.size());

    forAll(objectiveNames, objectivei)
    {
        const word& objectiveName = objectiveNames[objectivei];

        objectives_.set
        (
            objectivei,
            objective::New
            (
                mesh_,
                objectiveNamesDict.subDict(objectiveName),
                objectiveType,
                adjointSolverName,
                primalSolverName
            )
        );
    }

    if (objectives_.empty())
    {
        FatalIOErrorInFunction(objectiveNamesDict)
            << "No objectives have been set - cannot perform an optimisation"
            << exit(FatalIOError);
    }

    if (Pstream::master())
    {
        if (objectives_.size() > 1)
        {
            const Time& time = mesh_.time();
            weigthedObjectiveFile_.reset
            (
                new OFstream
                (
                    time.globalPath()/"optimisation"/"objective"
                        /time.timeName()/"weightedObjective"+adjointSolverName_
                )
            );

            unsigned int width = IOstream::defaultPrecision() + 5;
            weigthedObjectiveFile_()
                << setw(4) << "#" << " "
                << setw(width) << "weightedObjective" << " ";
            for (objective& objI : objectives_)
            {
                weigthedObjectiveFile_()
                    << setw(width) << objI.objectiveName() << " ";
            }
            weigthedObjectiveFile_()
                << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<objectiveManager> objectiveManager::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
{
    // Determine type of objectiveManager from objectiveType
    const word objectiveType(dict.get<word>("type"));
    const word managerType("objectiveManager" & objectiveType);

    auto* ctorPtr = dictionaryConstructorTable(managerType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "objectiveManagerType",
            managerType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<objectiveManager>
    (
        ctorPtr(mesh, dict, adjointSolverName, primalSolverName)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool objectiveManager::readDict(const dictionary& dict)
{
    for (objective& obj : objectives_)
    {
        obj.readDict
        (
            dict.subDict("objectiveNames").subDict(obj.objectiveName())
        );
    }

    return true;
}


void objectiveManager::updateNormalizationFactor()
{
    // Update normalization factors for all objectives
    for (objective& obj : objectives_)
    {
        if (obj.normalize())
        {
            obj.updateNormalizationFactor();
        }
    }
}


void objectiveManager::update()
{
    // Update all fields related to the objective function
    for (objective& obj : objectives_)
    {
        obj.update();
    }
}


void objectiveManager::updateOrNullify()
{
    // Update contributions to adjoint if true, otherwise return nulls
    for (objective& obj : objectives_)
    {
        if (obj.isWithinIntegrationTime())
        {
            obj.update();
        }
        else
        {
            obj.nullify();
        }
    }
}


void objectiveManager::incrementIntegrationTimes(const scalar timeSpan)
{
    // Update start and end integration times by adding the timeSpan
    // of the optimisation cycle
    for (objective& obj : objectives_)
    {
        obj.incrementIntegrationTimes(timeSpan);
    }
}


scalar objectiveManager::print()
{
    scalar objValue(Zero);
    Info<< "Adjoint solver " << adjointSolverName_ << endl;
    for (objective& obj : objectives_)
    {
        scalar cost = obj.JCycle();
        scalar weight = obj.weight();
        objValue += weight*cost;

        Info<< obj.objectiveName() << " : " << cost << endl;
    }

    Info<< "Weighted objective : " << objValue << nl << endl;

    return objValue;
}


bool objectiveManager::writeObjectives
(
    const scalar weightedObjective,
    const bool valid
)
{
    for (const objective& obj : objectives_)
    {
        // Write objective function to file
        obj.write();
        obj.writeMeanValue();
    }

    if (weigthedObjectiveFile_.valid())
    {
        unsigned int width = IOstream::defaultPrecision() + 5;
        weigthedObjectiveFile_()
            << setw(4) << mesh_.time().timeName() << " "
            << setw(width) << weightedObjective << " ";

        for (objective& objI : objectives_)
        {
            weigthedObjectiveFile_()
                << setw(width) << objI.JCycle() << " ";
        }
        weigthedObjectiveFile_() << endl;
    }

    return true;
}


void objectiveManager::updateAndWrite()
{
    updateNormalizationFactor();
    update();
    scalar weightedObjective = print();
    writeObjectives(weightedObjective);
}


PtrList<objective>& objectiveManager::getObjectiveFunctions()
{
    return objectives_;
}


const PtrList<objective>& objectiveManager::getObjectiveFunctions() const
{
    return objectives_;
}


const word& objectiveManager::adjointSolverName() const
{
    return adjointSolverName_;
}


const word& objectiveManager::primalSolverName() const
{
    return primalSolverName_;
}


void objectiveManager::checkIntegrationTimes() const
{
    for (const objective& obj : objectives_)
    {
        if (!obj.hasIntegrationStartTime() || !obj.hasIntegrationEndTime())
        {
            FatalErrorInFunction()
                << "Objective function " << obj.objectiveName()
                << " does not have a defined integration start or end time "
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
