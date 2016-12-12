/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "regionFunctionObject.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ObjectType>
bool Foam::functionObjects::regionFunctionObject::foundObject
(
    const word& fieldName
) const
{
    return obr().foundObject<ObjectType>(fieldName);
}


template<class ObjectType>
const ObjectType& Foam::functionObjects::regionFunctionObject::lookupObject
(
    const word& fieldName
) const
{
    return obr().lookupObject<ObjectType>(fieldName);
}


template<class ObjectType>
ObjectType& Foam::functionObjects::regionFunctionObject::lookupObjectRef
(
    const word& fieldName
) const
{
    return obr().lookupObjectRef<ObjectType>(fieldName);
}


template<class ObjectType>
const ObjectType* Foam::functionObjects::regionFunctionObject::lookupObjectPtr
(
    const word& fieldName
) const
{
    return obr().lookupObjectPtr<ObjectType>(fieldName);
}


template<class ObjectType>
ObjectType* Foam::functionObjects::regionFunctionObject::lookupObjectRefPtr
(
    const word& fieldName
) const
{
    return obr().lookupObjectRefPtr<ObjectType>(fieldName);
}


template<class ObjectType>
bool Foam::functionObjects::regionFunctionObject::store
(
    word& fieldName,
    const tmp<ObjectType>& tfield,
    bool cacheable
)
{
    if (cacheable && fieldName == tfield().name())
    {
        WarningInFunction
            << "Cannot store cache-able field with the name used in the cache."
            << nl
            << "    Either choose a different name or cache the field"
            << "    and use the 'writeObjects' functionObject."
            << endl;

        return false;
    }

    if (fieldName.size() && foundObject<ObjectType>(fieldName))
    {
        const ObjectType& field = lookupObject<ObjectType>(fieldName);

        // If there is a result field already registered, assign to the new
        // result field. Otherwise transfer ownership of the new result field to
        // the object registry
        if (&field != &tfield())
        {
            const_cast<ObjectType&>(field) = tfield;
        }
        else
        {
            obr().objectRegistry::store(tfield.ptr());
        }
    }
    else
    {
        if (fieldName.size() && fieldName != tfield().name())
        {
            tfield.ref().rename(fieldName);
        }
        else
        {
            fieldName = tfield().name();
        }

        obr().objectRegistry::store(tfield.ptr());
    }

    return true;
}


// ************************************************************************* //
