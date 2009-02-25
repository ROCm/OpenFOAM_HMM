/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Constructors & destructor for regIOobject.

\*---------------------------------------------------------------------------*/

#include "regIOobject.H"
#include "Time.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::regIOobject, 0);

int Foam::regIOobject::fileModificationSkew
(
    Foam::debug::optimisationSwitch("fileModificationSkew", 30)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from IOobject
Foam::regIOobject::regIOobject(const IOobject& io)
:
    IOobject(io),
    registered_(false),
    ownedByRegistry_(false),
    lastModified_(0),
    isPtr_(NULL)
{
    // Register with objectRegistry if requested
    if (registerObject())
    {
        checkIn();
    }
}


// Construct as copy
Foam::regIOobject::regIOobject(const regIOobject& rio)
:
    IOobject(rio),
    registered_(false),
    ownedByRegistry_(false),
    lastModified_(rio.lastModified_),
    isPtr_(NULL)
{
    // Do not register copy with objectRegistry
}


// Construct as copy, and transfering objectRegistry registration to copy
// if registerCopy is true
Foam::regIOobject::regIOobject(const regIOobject& rio, bool registerCopy)
:
    IOobject(rio),
    registered_(false),
    ownedByRegistry_(false),
    lastModified_(rio.lastModified_),
    isPtr_(NULL)
{
    if (registerCopy && rio.registered_)
    {
        const_cast<regIOobject&>(rio).checkOut();
        checkIn();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// Delete read stream, checkout from objectRegistry and destroy
Foam::regIOobject::~regIOobject()
{
    if (objectRegistry::debug)
    {
        Info<< "Destroying regIOobject called " << name()
            << " of type " << type()
            << " in directory " << path()
            << endl;
    }

    if (isPtr_)
    {
        delete isPtr_;
        isPtr_ = NULL;
    }

    // Check out of objectRegistry if not owned by the registry

    if (!ownedByRegistry_)
    {
        checkOut();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::regIOobject::checkIn()
{
    if (!registered_)
    {
        // multiple checkin of same object is disallowed - this would mess up
        // any mapping
        registered_ = db().checkIn(*this);

        // checkin on defaultRegion is allowed to fail, since subsetted meshes
        // are created with the same name as their originating mesh
        if (!registered_ && debug && name() != polyMesh::defaultRegion)
        {
            WarningIn("regIOobject::checkIn()")
                << "failed to register object " << objectPath()
                    << " the name already exists in the objectRegistry"
                << endl;
        }
    }

    return registered_;
}


bool Foam::regIOobject::checkOut()
{
    if (registered_)
    {
        registered_ = false;
        return db().checkOut(*this);
    }

    return false;
}


// Rename object and re-register with objectRegistry under new name
void Foam::regIOobject::rename(const word& newName)
{
    // TODO: verify that this works properly with unregistered objects
    // I suspect that it incorrectly registers them

    // Check out of objectRegistry
    checkOut();

    IOobject::rename(newName);

    // Re-register object with objectRegistry
    checkIn();
}


// Assign to IOobject
void Foam::regIOobject::operator=(const IOobject& io)
{
    if (isPtr_)
    {
        delete isPtr_;
        isPtr_ = NULL;
    }

    // Check out of objectRegistry
    checkOut();

    IOobject::operator=(io);

    // Re-register object with objectRegistry
    checkIn();
}


// ************************************************************************* //
