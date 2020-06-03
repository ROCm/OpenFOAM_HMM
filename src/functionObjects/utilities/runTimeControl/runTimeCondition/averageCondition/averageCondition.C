/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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

#include "averageCondition.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimeControls
{
    defineTypeNameAndDebug(averageCondition, 0);
    addToRunTimeSelectionTable(runTimeCondition, averageCondition, dictionary);
}
}
}

const Foam::Enum
<
    Foam::functionObjects::runTimeControls::averageCondition::windowType
>
Foam::functionObjects::runTimeControls::averageCondition::windowTypeNames
({
    { windowType::NONE, "none" },
    { windowType::APPROXIMATE, "approximate" },
    { windowType::EXACT, "exact" }
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimeControls::averageCondition::averageCondition
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    stateFunctionObject& state
)
:
    runTimeCondition(name, obr, dict, state),
    functionObjectName_(dict.get<word>("functionObject")),
    fieldNames_(dict.get<wordList>("fields")),
    tolerance_(dict.get<scalar>("tolerance")),
    window_(dict.getOrDefault<scalar>("window", -1)),
    windowType_
    (
        window_ > 0
      ? windowTypeNames.get("windowType", dict)
      : windowType::NONE
    ),
    totalTime_(fieldNames_.size(), scalar(0)),
    resetOnRestart_(dict.getOrDefault("resetOnRestart", false)),
    nIterStartUp_(dict.getOrDefault<label>("nIterStartUp", 10)),
    iter_(-1)
{
    dictionary& conditionDict = this->conditionDict();

    forAll(fieldNames_, fieldi)
    {
        const word& fieldName = fieldNames_[fieldi];

        if (resetOnRestart_)
        {
            conditionDict.set(fieldName, dictionary());
        }
        else
        {
            if (conditionDict.found(fieldName))
            {
                const dictionary& valueDict = conditionDict.subDict(fieldName);
                valueDict.readIfPresent("totalTime", totalTime_[fieldi]);
            }
            else
            {
                conditionDict.set(fieldName, dictionary());
            }
        }
    }

    conditionDict.readIfPresent("iter", iter_);
}


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

bool Foam::functionObjects::runTimeControls::averageCondition::apply()
{
    if (!active_)
    {
        return true;
    }

    bool satisfied = iter_ > nIterStartUp_;

    ++iter_;

    const scalar dt = obr_.time().deltaTValue();

    Log << "    " << type() << ": " << name_ << " averages:" << nl;

    DynamicList<label> unprocessedFields(fieldNames_.size());

    forAll(fieldNames_, fieldi)
    {
        totalTime_[fieldi] += dt;

        bool processed = false;
        calc<scalar>(fieldi, satisfied, processed);
        calc<vector>(fieldi, satisfied, processed);
        calc<sphericalTensor>(fieldi, satisfied, processed);
        calc<symmTensor>(fieldi, satisfied, processed);
        calc<tensor>(fieldi, satisfied, processed);

        if (!processed)
        {
            unprocessedFields.append(fieldi);
        }
    }

    if (unprocessedFields.size())
    {
        WarningInFunction
            << "From function object: " << functionObjectName_ << nl
            << "Unprocessed fields:" << nl;

        for (const label fieldi : unprocessedFields)
        {
            Info<< "        " << fieldNames_[fieldi] << nl;
        }

        if (unprocessedFields.size() == fieldNames_.size())
        {
            satisfied = false;
        }
    }

    Log << endl;

    return satisfied;
}


void Foam::functionObjects::runTimeControls::averageCondition::write()
{
    dictionary& conditionDict = this->conditionDict();

    forAll(fieldNames_, fieldi)
    {
        const word& fieldName = fieldNames_[fieldi];

        if (conditionDict.found(fieldName))
        {
            dictionary& valueDict = conditionDict.subDict(fieldName);
            valueDict.add("totalTime", totalTime_[fieldi], true);
        }
        else
        {
            dictionary valueDict;
            valueDict.add("totalTime", totalTime_[fieldi], true);
            conditionDict.add(fieldName, valueDict);
        }
    }

    conditionDict.set("iter", iter_);
}


// ************************************************************************* //
