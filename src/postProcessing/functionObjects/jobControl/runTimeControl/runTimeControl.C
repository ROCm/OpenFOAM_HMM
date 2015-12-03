/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(runTimeControl, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::runTimeControl::runTimeControl
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjectState(obr, name),
    obr_(obr),
    conditions_(),
    groupMap_(),
    nWriteStep_(0),
    writeStepI_(0)
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (setActive<fvMesh>())
    {
        read(dict);

        // Check that some conditions are set
        if (conditions_.empty())
        {
            Info<< type() << " " << name_ << " output:" << nl
                << "    No conditions present - deactivating" << nl
                << endl;

            active_ = false;
        }
        else
        {
            // Check that at least one condition is active
            active_ = false;
            forAll(conditions_, conditionI)
            {
                if (conditions_[conditionI].active())
                {
                    active_ = true;
                    break;
                }
            }

            if (!active_)
            {
                Info<< type() << " " << name_ << " output:" << nl
                    << "    All conditions inactive - deactivating" << nl
                    << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::runTimeControl::~runTimeControl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::runTimeControl::read(const dictionary& dict)
{
    if (active_)
    {
        const dictionary& conditionsDict = dict.subDict("conditions");
        const wordList conditionNames(conditionsDict.toc());
        conditions_.setSize(conditionNames.size());

        label uniqueGroupI = 0;
        forAll(conditionNames, conditionI)
        {
            const word& conditionName = conditionNames[conditionI];
            const dictionary& dict = conditionsDict.subDict(conditionName);

            conditions_.set
            (
                conditionI,
                runTimeCondition::New(conditionName, obr_, dict, *this)
            );

            label groupI = conditions_[conditionI].groupID();

            if (groupMap_.insert(groupI, uniqueGroupI))
            {
                uniqueGroupI++;
            }
        }

        dict.readIfPresent("nWriteStep", nWriteStep_);
    }
}


void Foam::runTimeControl::execute()
{
    if (!active_)
    {
        return;
    }

    Info<< type() << " " << name_ << " output:" << nl;

    // IDs of satisfied conditions
    DynamicList<label> IDs(conditions_.size());

    // Run stops only if all conditions within a group are satisfied
    List<bool> groupSatisfied(groupMap_.size(), true);
    List<bool> groupActive(groupMap_.size(), false);

    forAll(conditions_, conditionI)
    {
        runTimeCondition& condition = conditions_[conditionI];

        if (condition.active())
        {
            bool conditionSatisfied = condition.apply();

            label groupI = condition.groupID();

            Map<label>::const_iterator conditionIter = groupMap_.find(groupI);

            if (conditionIter == groupMap_.end())
            {
                FatalErrorIn("void Foam::runTimeControl::execute()")
                    << "group " << groupI << " not found in map"
                    << abort(FatalError);
            }

            if (conditionSatisfied)
            {
                IDs.append(conditionI);

                groupActive[conditionIter()] = true;

                if (groupI == -1)
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
    forAll(groupSatisfied, groupI)
    {
        if (groupSatisfied[groupI] && groupActive[groupI])
        {
            done = true;
            break;
        }
    }

    if (done)
    {
        forAll(IDs, conditionI)
        {
            Info<< "    " << conditions_[conditionI].type() << ": "
                <<  conditions_[conditionI].name()
                << " condition satisfied" << nl;
        }


        // Set to write a data dump or finalise the calculation
        Time& time = const_cast<Time&>(obr_.time());

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
}


void Foam::runTimeControl::end()
{
    // Do nothing
}


void Foam::runTimeControl::timeSet()
{
    // Do nothing
}


void Foam::runTimeControl::write()
{
    if (active_)
    {
        forAll(conditions_, conditionI)
        {
            conditions_[conditionI].write();
        }
    }
}


// ************************************************************************* //
