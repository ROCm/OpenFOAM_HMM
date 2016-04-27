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

#include "functionObjectState.H"
#include "Time.H"

const Foam::word Foam::functionObjectState::resultsName_ = "results";

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const Foam::IOdictionary& Foam::functionObjectState::stateDict() const
{
    return obr_.time().functionObjects().stateDict();
}


Foam::IOdictionary& Foam::functionObjectState::stateDict()
{
    return const_cast<IOdictionary&>(obr_.time().functionObjects().stateDict());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjectState::functionObjectState
(
    const objectRegistry& obr,
    const word& name
)
:
    obr_(obr),
    name_(name),
    active_(true)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjectState::~functionObjectState()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::functionObjectState::name() const
{
    return name_;
}


bool Foam::functionObjectState::active() const
{
    return active_;
}


Foam::dictionary& Foam::functionObjectState::propertyDict()
{
    IOdictionary& stateDict = this->stateDict();

    if (!stateDict.found(name_))
    {
        stateDict.add(name_, dictionary());
    }

    return stateDict.subDict(name_);
}


bool Foam::functionObjectState::foundProperty(const word& entryName) const
{
    const IOdictionary& stateDict = this->stateDict();

    if (stateDict.found(name_))
    {
        const dictionary& baseDict = stateDict.subDict(name_);
        return baseDict.found(entryName);
    }

    return false;
}


Foam::word Foam::functionObjectState::resultType(const word& entryName) const
{
    return objectResultType(name_, entryName);
}


Foam::word Foam::functionObjectState::objectResultType
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

            forAllConstIter(dictionary, objectDict, iter)
            {
                const dictionary& dict = iter().dict();

                if (dict.found(entryName))
                {
                    return dict.dictName();
                }
            }
        }
    }

    return result;
}


Foam::List<Foam::word> Foam::functionObjectState::objectResultEntries() const
{
    return objectResultEntries(name_);
}


Foam::List<Foam::word> Foam::functionObjectState::objectResultEntries
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

            forAllConstIter(dictionary, objectDict, iter)
            {
                const dictionary& dict = iter().dict();
                result.append(dict.toc());
            }
        }
    }

    wordList entries;
    entries.transfer(result);

    return entries;
}


// ************************************************************************* //
