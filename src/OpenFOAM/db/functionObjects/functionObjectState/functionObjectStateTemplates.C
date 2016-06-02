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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjectState::setActive()
{
    active_ = true;

    if (!isA<Type>(obr_))
    {
        WarningInFunction
            << "No " << Type::typeName << " available, deactivating " << name_
            << endl;

        active_ = false;
    }

    return active_;
}


template<class Type>
Type Foam::functionObjectState::getProperty
(
    const word& entryName,
    const Type& defaultValue
) const
{
    Type result = defaultValue;
    getProperty(entryName, result);
    return result;
}


template<class Type>
void Foam::functionObjectState::getProperty
(
    const word& entryName,
    Type& value
) const
{
    getObjectProperty(name_, entryName, value);
}


template<class Type>
void Foam::functionObjectState::setProperty
(
    const word& entryName,
    const Type& value
)
{
    setObjectProperty(name_, entryName, value);
}


template<class Type>
Type Foam::functionObjectState::getObjectProperty
(
    const word& objectName,
    const word& entryName,
    const Type& defaultValue
) const
{
    Type result = defaultValue;
    getObjectProperty(objectName, entryName, result);
    return result;
}


template<class Type>
void Foam::functionObjectState::getObjectProperty
(
    const word& objectName,
    const word& entryName,
    Type& value
) const
{
    const IOdictionary& stateDict = this->stateDict();

    if (stateDict.found(objectName))
    {
        const dictionary& baseDict = stateDict.subDict(objectName);
        if (baseDict.found(entryName))
        {
            if (baseDict.isDict(entryName))
            {
                value = baseDict.subDict(entryName);
            }
            else
            {
                baseDict.lookup(entryName) >> value;
            }
        }
    }
}


template<class Type>
void Foam::functionObjectState::setObjectProperty
(
    const word& objectName,
    const word& entryName,
    const Type& value
)
{
    IOdictionary& stateDict = this->stateDict();

    if (!stateDict.found(objectName))
    {
        stateDict.add(objectName, dictionary());
    }

    dictionary& baseDict = stateDict.subDict(objectName);
    baseDict.add(entryName, value, true);
}


template<class Type>
void Foam::functionObjectState::setResult
(
    const word& entryName,
    const Type& value
)
{
    setObjectResult(name_, entryName, value);
}


template<class Type>
void Foam::functionObjectState::setObjectResult
(
    const word& objectName,
    const word& entryName,
    const Type& value
)
{
    IOdictionary& stateDict = this->stateDict();

    if (!stateDict.found(resultsName_))
    {
        stateDict.add(resultsName_, dictionary());
    }

    dictionary& resultsDict = stateDict.subDict(resultsName_);

    if (!resultsDict.found(objectName))
    {
        resultsDict.add(name_, dictionary());
    }

    dictionary& objectDict = resultsDict.subDict(objectName);

    const word& dictTypeName = pTraits<Type>::typeName;

    if (!objectDict.found(dictTypeName))
    {
        objectDict.add(dictTypeName, dictionary());
    }

    dictionary& resultTypeDict = objectDict.subDict(dictTypeName);

    resultTypeDict.add(entryName, value, true);
}


template<class Type>
Type Foam::functionObjectState::getResult
(
    const word& entryName,
    const Type& defaultValue
) const
{
    return getObjectResult(name_, entryName, defaultValue);
}


template<class Type>
Type Foam::functionObjectState::getObjectResult
(
    const word& objectName,
    const word& entryName,
    const Type& defaultValue
) const
{
    Type result = defaultValue;
    getObjectResult(objectName, entryName, result);
    return result;
}


template<class Type>
void Foam::functionObjectState::getObjectResult
(
    const word& objectName,
    const word& entryName,
    Type& value
) const
{
    const IOdictionary& stateDict = this->stateDict();

    if (stateDict.found(resultsName_))
    {
        const dictionary& resultsDict = stateDict.subDict(resultsName_);

        if (resultsDict.found(objectName))
        {
            const dictionary& objectDict = resultsDict.subDict(objectName);

            const word& dictTypeName = pTraits<Type>::typeName;

            if (objectDict.found(dictTypeName))
            {
                const dictionary& resultTypeDict =
                    objectDict.subDict(dictTypeName);

                resultTypeDict.readIfPresent<Type>(entryName, value);
            }
        }
    }
}


// ************************************************************************* //
