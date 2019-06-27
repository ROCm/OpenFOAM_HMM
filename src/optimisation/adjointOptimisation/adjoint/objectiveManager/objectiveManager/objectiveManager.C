/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2007-2019 PCOpt/NTUA
                            | Copyright (C) 2013-2019 FOSS GP
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
    objectives_(0)
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

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(managerType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown objectiveManagerType type " << managerType
            << nl << nl
            << "Valid objectiveManagerTypes are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<objectiveManager>
    (
        cstrIter()(mesh, dict, adjointSolverName, primalSolverName)
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
        obj.updateNormalizationFactor();
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


scalar objectiveManager::print()
{
    scalar objValue(Zero);
    for (objective& obj : objectives_)
    {
        scalar cost = obj.J();
        scalar weight = obj.weight();
        objValue += weight*cost;

        Info<< obj.type() << " : " << cost << endl;
    }

    Info<< "Objective function manager" << nl
        << "    Weighted Lagrangian " << " : " << objValue << nl << endl;

    return objValue;
}


bool objectiveManager::write(const bool valid) const
{
    for (const objective& obj : objectives_)
    {
        // Write objective function to file
        obj.write();
        obj.writeMeanValue();
    }

    return true;
}


void objectiveManager::updateAndWrite()
{
    updateNormalizationFactor();
    update();
    print();
    write();
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
