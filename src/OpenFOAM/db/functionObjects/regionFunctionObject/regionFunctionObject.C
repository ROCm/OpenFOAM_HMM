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
#include "Time.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(regionFunctionObject, 0);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

const Foam::objectRegistry&
Foam::functionObjects::regionFunctionObject::whichSubRegistry
(
    const objectRegistry& obr,
    const dictionary& dict
)
{
    word subName;
    if (dict.readIfPresent("subRegion", subName))
    {
        return obr.lookupObject<objectRegistry>(subName);
    }
    else
    {
        return obr;
    }
}


const Foam::objectRegistry&
Foam::functionObjects::regionFunctionObject::obr() const
{
    return subObr_;
}


bool Foam::functionObjects::regionFunctionObject::writeObject
(
    const word& fieldName
)
{
    const regIOobject* objPtr =
        this->lookupObjectPtr<regIOobject>(fieldName);

    if (objPtr)
    {
        Log << "    functionObjects::" << type() << " " << name()
            << " writing field: " << objPtr->name() << endl;

        objPtr->write();

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::functionObjects::regionFunctionObject::clearObject
(
    const word& fieldName
)
{
    regIOobject* objPtr = lookupObjectRefPtr<regIOobject>(fieldName);
    if (objPtr)
    {
        if (objPtr->ownedByRegistry())
        {
            return objPtr->checkOut();
        }
        else
        {
            return false;
        }
    }
    else
    {
        return true;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::regionFunctionObject::regionFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    stateFunctionObject(name, runTime),
    obr_
    (
        runTime.lookupObject<objectRegistry>
        (
            dict.lookupOrDefault("region", polyMesh::defaultRegion)
        )
    ),
    subObr_(whichSubRegistry(obr_, dict))
{}


Foam::functionObjects::regionFunctionObject::regionFunctionObject
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    stateFunctionObject(name, obr.time()),
    obr_(obr),
    subObr_(whichSubRegistry(obr_, dict))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::regionFunctionObject::~regionFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::regionFunctionObject::read(const dictionary& dict)
{
    return stateFunctionObject::read(dict);
}


// ************************************************************************* //
