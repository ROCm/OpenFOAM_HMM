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

#include "stateFunctionObject.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::functionObjects::stateFunctionObject::resultsName_ =
    "results";


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const Foam::IOdictionary&
Foam::functionObjects::stateFunctionObject::stateDict() const
{
    return time_.functionObjects().stateDict();
}


Foam::IOdictionary& Foam::functionObjects::stateFunctionObject::stateDict()
{
    return const_cast<IOdictionary&>(time_.functionObjects().stateDict());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::stateFunctionObject::stateFunctionObject
(
    const word& name,
    const Time& runTime
)
:
    timeFunctionObject(name, runTime)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dictionary& Foam::functionObjects::stateFunctionObject::propertyDict()
{
    IOdictionary& stateDict = this->stateDict();

    if (!stateDict.found(name()))
    {
        stateDict.add(name(), dictionary());
    }

    return stateDict.subDict(name());
}


bool Foam::functionObjects::stateFunctionObject::setTrigger
(
    const label triggeri
)
{
    IOdictionary& stateDict = this->stateDict();

    label oldTriggeri =
        stateDict.getOrDefault<label>("triggerIndex", labelMin);

    if (triggeri > oldTriggeri)
    {
        stateDict.set("triggerIndex", triggeri);
        return true;
    }

    return false;
}


Foam::label Foam::functionObjects::stateFunctionObject::getTrigger() const
{
    const IOdictionary& stateDict = this->stateDict();

    return stateDict.getOrDefault<label>("triggerIndex", labelMin);
}


bool Foam::functionObjects::stateFunctionObject::foundProperty
(
    const word& entryName
) const
{
    const IOdictionary& stateDict = this->stateDict();

    if (stateDict.found(name()))
    {
        const dictionary& baseDict = stateDict.subDict(name());
        return baseDict.found(entryName);
    }

    return false;
}


bool Foam::functionObjects::stateFunctionObject::getDict
(
    const word& entryName,
    dictionary& dict
) const
{
    return getObjectDict(name(), entryName, dict);
}


bool Foam::functionObjects::stateFunctionObject::getObjectDict
(
    const word& objectName,
    const word& entryName,
    dictionary& dict
) const
{
    const IOdictionary& stateDict = this->stateDict();

    if (stateDict.found(objectName))
    {
        const dictionary& baseDict = stateDict.subDict(objectName);
        if (baseDict.found(entryName) && baseDict.isDict(entryName))
        {
            dict = baseDict.subDict(entryName);
            return true;
        }
    }

    return false;
}


Foam::word Foam::functionObjects::stateFunctionObject::resultType
(
    const word& entryName
) const
{
    return objectResultType(name(), entryName);
}


Foam::word Foam::functionObjects::stateFunctionObject::objectResultType
(
    const word& objectName,
    const word& entryName
) const
{
    word result = word::null;
    const IOdictionary& stateDict = this->stateDict();

    if (stateDict.found(resultsName_))
    {
        const dictionary& resultsDict = stateDict.subDict(resultsName_);

        if (resultsDict.found(objectName))
        {
            const dictionary& objectDict = resultsDict.subDict(objectName);

            for (const entry& dEntry : objectDict)
            {
                const dictionary& dict = dEntry.dict();

                if (dict.found(entryName))
                {
                    return dict.dictName();
                }
            }
        }
    }

    return result;
}


Foam::List<Foam::word>
Foam::functionObjects::stateFunctionObject::objectResultEntries() const
{
    return objectResultEntries(name());
}


Foam::List<Foam::word>
Foam::functionObjects::stateFunctionObject::objectResultEntries
(
    const word& objectName
) const
{
    DynamicList<word> result(2);

    const IOdictionary& stateDict = this->stateDict();

    if (stateDict.found(resultsName_))
    {
        const dictionary& resultsDict = stateDict.subDict(resultsName_);

        if (resultsDict.found(objectName))
        {
            const dictionary& objectDict = resultsDict.subDict(objectName);

            for (const entry& dEntry : objectDict)
            {
                const dictionary& dict = dEntry.dict();

                result.append(dict.toc());
            }
        }
    }

    wordList entries;
    entries.transfer(result);

    return entries;
}

void Foam::functionObjects::stateFunctionObject::writeResultEntries
(
    Ostream& os
) const
{
    writeResultEntries(name(), os);
}


void Foam::functionObjects::stateFunctionObject::writeResultEntries
(
    const word& objectName,
    Ostream& os
) const
{
    const IOdictionary& stateDict = this->stateDict();

    if (stateDict.found(resultsName_))
    {
        const dictionary& resultsDict = stateDict.subDict(resultsName_);

        if (resultsDict.found(objectName))
        {
            const dictionary& objectDict = resultsDict.subDict(objectName);

            for (const word& dataFormat : objectDict.sortedToc())
            {
                os  << "    Type: " << dataFormat << nl;

                const dictionary& resultDict = objectDict.subDict(dataFormat);

                for (const word& result : resultDict.sortedToc())
                {
                    os << "        " << result << nl;
                }
            }
        }
    }
}


void Foam::functionObjects::stateFunctionObject::writeAllResultEntries
(
    Ostream& os
) const
{
    const IOdictionary& stateDict = this->stateDict();

    if (stateDict.found(resultsName_))
    {
        const dictionary& resultsDict = stateDict.subDict(resultsName_);

        const wordList allObjectNames = resultsDict.sortedToc();

        for (const word& objectName : allObjectNames)
        {
            os  << "Object: " << objectName << endl;

            writeResultEntries(objectName, os);
        }
    }
}


// ************************************************************************* //
