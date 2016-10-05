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

#include "runTimeControl.H"
#include "dictionary.H"
#include "runTimeCondition.H"
#include "fvMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimeControls
{
    defineTypeNameAndDebug(runTimeControl, 0);
    addToRunTimeSelectionTable(functionObject, runTimeControl, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimeControls::runTimeControl::runTimeControl
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    conditions_(),
    groupMap_(),
    nWriteStep_(0),
    writeStepI_(0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::runTimeControls::runTimeControl::~runTimeControl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::runTimeControls::runTimeControl::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);

    const dictionary& conditionsDict = dict.subDict("conditions");
    const wordList conditionNames(conditionsDict.toc());
    conditions_.setSize(conditionNames.size());

    label uniqueGroupi = 0;
    forAll(conditionNames, conditioni)
    {
        const word& conditionName = conditionNames[conditioni];
        const dictionary& dict = conditionsDict.subDict(conditionName);

        conditions_.set
        (
            conditioni,
            runTimeCondition::New(conditionName, obr_, dict, *this)
        );

        label groupi = conditions_[conditioni].groupID();

        if (groupMap_.insert(groupi, uniqueGroupi))
        {
            uniqueGroupi++;
        }
    }

    dict.readIfPresent("nWriteStep", nWriteStep_);

    // Check that some conditions are set
    if (conditions_.empty())
    {
        Info<< type() << " " << name() << " output:" << nl
            << "    No conditions present" << nl
            << endl;
    }
    else
    {
        // Check that at least one condition is active
        bool active = false;
        forAll(conditions_, conditioni)
        {
            if (conditions_[conditioni].active())
            {
                active = true;
                break;
            }
        }

        if (!active)
        {
            Info<< type() << " " << name() << " output:" << nl
                << "    All conditions are inactive" << nl
                << endl;
        }
    }

    return true;
}


bool Foam::functionObjects::runTimeControls::runTimeControl::execute()
{
    Info<< type() << " " << name() << " output:" << nl;

    // IDs of satisfied conditions
    DynamicList<label> IDs(conditions_.size());

    // Run stops only if all conditions within a group are satisfied
    List<bool> groupSatisfied(groupMap_.size(), true);
    List<bool> groupActive(groupMap_.size(), false);

    forAll(conditions_, conditioni)
    {
        runTimeCondition& condition = conditions_[conditioni];

        if (condition.active())
        {
            bool conditionSatisfied = condition.apply();

            label groupi = condition.groupID();

            Map<label>::const_iterator conditionIter = groupMap_.find(groupi);

            if (conditionIter == groupMap_.end())
            {
                FatalErrorInFunction
                    << "group " << groupi << " not found in map"
                    << abort(FatalError);
            }

            if (conditionSatisfied)
            {
                IDs.append(conditioni);

                groupActive[conditionIter()] = true;

                if (groupi == -1)
                {
                    // Condition not part of a group - only requires this to be
                    // satisfied for completion flag to be set
                    groupSatisfied[conditionIter()] = true;
                    break;
                }
            }
            else
            {
                groupSatisfied[conditionIter()] = false;
            }
        }
    }

    bool done = false;
    forAll(groupSatisfied, groupi)
    {
        if (groupSatisfied[groupi] && groupActive[groupi])
        {
            done = true;
            break;
        }
    }

    if (done)
    {
        forAll(IDs, conditioni)
        {
            Info<< "    " << conditions_[conditioni].type() << ": "
                <<  conditions_[conditioni].name()
                << " condition satisfied" << nl;
        }


        // Set to write a data dump or finalise the calculation
        Time& time = const_cast<Time&>(time_);

        if (writeStepI_ < nWriteStep_ - 1)
        {
            writeStepI_++;
            Info<< "    Writing fields - step " << writeStepI_ << nl;
            time.writeNow();
        }
        else
        {
            Info<< "    Stopping calculation" << nl
                << "    Writing fields - final step" << nl;
            time.writeAndEnd();
        }
    }
    else
    {
        Info<< "    Conditions not met - calculations proceeding" << nl;
    }

    Info<< endl;

    return true;
}


bool Foam::functionObjects::runTimeControls::runTimeControl::write()
{
    forAll(conditions_, conditioni)
    {
        conditions_[conditioni].write();
    }

    return true;
}


// ************************************************************************* //
