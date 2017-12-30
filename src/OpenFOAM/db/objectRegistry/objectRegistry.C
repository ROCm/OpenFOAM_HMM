/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2017 OpenCFD Ltd.
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

#include "objectRegistry.H"
#include "Time.H"
#include "predicates.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(objectRegistry, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::objectRegistry::parentNotTime() const
{
    return (&parent_ != dynamic_cast<const objectRegistry*>(&time_));
}


// * * * * * * * * * * * * * * * * Constructors *  * * * * * * * * * * * * * //

Foam::objectRegistry::objectRegistry
(
    const Time& t,
    const label nIoObjects
)
:
    regIOobject
    (
        IOobject
        (
            word::validate(t.caseName()),
            t.path(),
            t,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        ),
        true    // to flag that this is the top-level regIOobject
    ),
    HashTable<regIOobject*>(nIoObjects),
    time_(t),
    parent_(t),
    dbDir_(name()),
    event_(1)
{}


Foam::objectRegistry::objectRegistry
(
    const IOobject& io,
    const label nIoObjects
)
:
    regIOobject(io),
    HashTable<regIOobject*>(nIoObjects),
    time_(io.time()),
    parent_(io.db()),
    dbDir_(parent_.dbDir()/local()/name()),
    event_(1)
{
    writeOpt() = IOobject::AUTO_WRITE;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::objectRegistry::~objectRegistry()
{
    List<regIOobject*> myObjects(size());
    label nMyObjects = 0;

    for (iterator iter = begin(); iter != end(); ++iter)
    {
        if (iter()->ownedByRegistry())
        {
            myObjects[nMyObjects++] = iter();
        }
    }

    for (label i=0; i < nMyObjects; ++i)
    {
        checkOut(*myObjects[i]);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::HashTable<Foam::wordHashSet> Foam::objectRegistry::classes() const
{
    return classesImpl(*this, predicates::always());
}


Foam::HashTable<Foam::wordHashSet>
Foam::objectRegistry::classes(const wordRe& matcher) const
{
    return classesImpl(*this, matcher);
}


Foam::HashTable<Foam::wordHashSet>
Foam::objectRegistry::classes(const wordRes& matcher) const
{
    return classesImpl(*this, matcher);
}


Foam::wordList Foam::objectRegistry::names() const
{
    return HashTable<regIOobject*>::toc();
}


Foam::wordList Foam::objectRegistry::sortedNames() const
{
    return HashTable<regIOobject*>::sortedToc();
}


Foam::wordList Foam::objectRegistry::names(const word& clsName) const
{
    wordList objNames(size());

    label count=0;
    for (const_iterator iter = cbegin(); iter != cend(); ++iter)
    {
        if (iter()->type() == clsName)
        {
            objNames[count++] = iter.key();
        }
    }

    objNames.setSize(count);

    return objNames;
}


Foam::wordList Foam::objectRegistry::sortedNames(const word& clsName) const
{
    wordList objNames = names(clsName);
    Foam::sort(objNames);

    return objNames;
}


const Foam::objectRegistry& Foam::objectRegistry::subRegistry
(
    const word& name,
    const bool forceCreate,
    const bool recursive
) const
{
    if (forceCreate && !foundObject<objectRegistry>(name, recursive))
    {
        objectRegistry* subObr = new objectRegistry
        (
            IOobject
            (
                name,
                time().constant(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );
        subObr->store();
    }

    return lookupObject<objectRegistry>(name, recursive);
}


Foam::label Foam::objectRegistry::getEvent() const
{
    label curEvent = event_++;

    if (event_ == labelMax)
    {
        if (objectRegistry::debug)
        {
            WarningInFunction
                << "Event counter has overflowed. "
                << "Resetting counter on all dependent objects." << nl
                << "This might cause extra evaluations." << endl;
        }

        // Reset event counter
        curEvent = 1;
        event_ = 2;


        // No need to reset dependent objects; overflow is now handled
        // in regIOobject::upToDate
    }

    return curEvent;
}


bool Foam::objectRegistry::checkIn(regIOobject& io) const
{
    if (objectRegistry::debug)
    {
        Pout<< "objectRegistry::checkIn(regIOobject&) : "
            << name() << " : checking in " << io.name()
            << " of type " << io.type()
            << endl;
    }

    return const_cast<objectRegistry&>(*this).insert(io.name(), &io);
}


bool Foam::objectRegistry::checkOut(regIOobject& io) const
{
    iterator iter = const_cast<objectRegistry&>(*this).find(io.name());

    if (iter.found())
    {
        if (objectRegistry::debug)
        {
            Pout<< "objectRegistry::checkOut(regIOobject&) : "
                << name() << " : checking out " << iter.key()
                << endl;
        }

        if (iter() != &io)
        {
            if (objectRegistry::debug)
            {
                WarningInFunction
                    << name() << " : attempt to checkOut copy of "
                    << iter.key()
                    << endl;
            }

            return false;
        }
        else
        {
            regIOobject* object = iter();

            bool hasErased = const_cast<objectRegistry&>(*this).erase(iter);

            if (io.ownedByRegistry())
            {
                delete object;
            }

            return hasErased;
        }
    }
    else
    {
        if (objectRegistry::debug)
        {
            Pout<< "objectRegistry::checkOut(regIOobject&) : "
                << name() << " : could not find " << io.name()
                << " in registry " << name()
                << endl;
        }
    }

    return false;
}


void Foam::objectRegistry::rename(const word& newName)
{
    regIOobject::rename(newName);

    // adjust dbDir_ as well
    string::size_type i = dbDir_.rfind('/');

    if (i == string::npos)
    {
        dbDir_ = newName;
    }
    else
    {
        dbDir_.replace(i+1, string::npos, newName);
    }
}


bool Foam::objectRegistry::modified() const
{
    forAllConstIter(HashTable<regIOobject*>, *this, iter)
    {
        if (iter()->modified())
        {
            return true;
        }
    }

    return false;
}


void Foam::objectRegistry::readModifiedObjects()
{
    for (iterator iter = begin(); iter != end(); ++iter)
    {
        if (objectRegistry::debug)
        {
            Pout<< "objectRegistry::readModifiedObjects() : "
                << name() << " : Considering reading object "
                << iter.key() << endl;
        }

        iter()->readIfModified();
    }
}


bool Foam::objectRegistry::readIfModified()
{
    readModifiedObjects();
    return true;
}


bool Foam::objectRegistry::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool valid
) const
{
    bool ok = true;

    forAllConstIter(HashTable<regIOobject*>, *this, iter)
    {
        if (objectRegistry::debug)
        {
            Pout<< "objectRegistry::write() : "
                << name() << " : Considering writing object "
                << iter.key()
                << " of type " << iter()->type()
                << " with writeOpt " << iter()->writeOpt()
                << " to file " << iter()->objectPath()
                << endl;
        }

        if (iter()->writeOpt() != NO_WRITE)
        {
            ok = iter()->writeObject(fmt, ver, cmp, valid) && ok;
        }
    }

    return ok;
}


// ************************************************************************* //
