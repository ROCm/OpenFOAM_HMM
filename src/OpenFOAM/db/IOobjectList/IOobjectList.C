/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "IOobjectList.H"
#include "Time.H"
#include "IOList.H"
#include "predicates.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::IOobjectList::checkNames(wordList& masterNames, const bool syncPar)
{
    // Sort for consistent order on all processors
    Foam::sort(masterNames);

    if (syncPar && Pstream::parRun())
    {
        const wordList localNames(masterNames);
        Pstream::scatter(masterNames);

        if (localNames != masterNames)
        {
            FatalErrorInFunction
                << "Objects not synchronised across processors." << nl
                << "Master has " << flatOutput(masterNames) << nl
                << "Processor " << Pstream::myProcNo()
                << " has " << flatOutput(localNames)
                << exit(FatalError);

            return false;
        }
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IOobjectList::IOobjectList()
:
    HashPtrTable<IOobject>()
{}


Foam::IOobjectList::IOobjectList(const label nObjects)
:
    HashPtrTable<IOobject>(nObjects)  // Could also use 2*nObjects instead
{}


Foam::IOobjectList::IOobjectList(const IOobjectList& list)
:
    HashPtrTable<IOobject>(list)
{}


Foam::IOobjectList::IOobjectList(IOobjectList&& list)
:
    HashPtrTable<IOobject>(std::move(list))
{}


Foam::IOobjectList::IOobjectList
(
    const objectRegistry& db,
    const fileName& instance,
    const fileName& local,
    IOobject::readOption r,
    IOobject::writeOption w,
    bool registerObject
)
:
    HashPtrTable<IOobject>()
{
    word newInstance;
    fileNameList objNames = fileHandler().readObjects
    (
        db,
        instance,
        local,
        newInstance
    );

    for (const auto& objName : objNames)
    {
        auto objectPtr = autoPtr<IOobject>::New
        (
            objName,
            newInstance,
            local,
            db,
            r,
            w,
            registerObject
        );

        bool ok = false;
        const bool throwingIOerr = FatalIOError.throwExceptions();

        try
        {
            // Use object with local scope and current instance (no searching)
            ok = objectPtr->typeHeaderOk<IOList<label>>(false, false);
        }
        catch (Foam::IOerror& err)
        {
            Warning
                << err << nl << endl;
        }

        FatalIOError.throwExceptions(throwingIOerr);

        if (ok)
        {
            insert(objectPtr->name(), objectPtr);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::IOobjectList::add(autoPtr<IOobject>& objectPtr)
{
    if (objectPtr.valid())
    {
        return insert(objectPtr->name(), objectPtr);
    }

    return false;
}


bool Foam::IOobjectList::add(autoPtr<IOobject>&& objectPtr)
{
    if (objectPtr.valid())
    {
        return insert(objectPtr->name(), objectPtr);
    }

    return false;
}


Foam::label Foam::IOobjectList::append(const IOobjectList& other)
{
    label count = 0;

    forAllConstIters(other, iter)
    {
        if (!found(iter.key()))
        {
            if (IOobject::debug)
            {
                InfoInFunction << "Copy append " << iter.key() << nl;
            }

            set(iter.key(), new IOobject(*(iter.object())));
            ++count;
        }
    }

    return count;
}


Foam::label Foam::IOobjectList::append(IOobjectList&& other)
{
    // Remove by name to avoid uncertainties about invalid iterators

    label count = 0;

    wordList keys(other.toc());

    for (const word& key : keys)
    {
        if (!found(key))
        {
            if (IOobject::debug)
            {
                InfoInFunction << "Move append " << key << nl;
            }

            if (add(other.remove(key)))
            {
                ++count;
            }
        }
    }

    return count;
}


bool Foam::IOobjectList::remove(const IOobject& io)
{
    return erase(io.name());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::IOobject* Foam::IOobjectList::findObject(const word& objName) const
{
    const_iterator iter = cfind(objName);

    if (iter.found())
    {
        if (IOobject::debug)
        {
            InfoInFunction << "Found " << objName << endl;
        }

        return const_cast<IOobject*>(*iter);
    }

    if (IOobject::debug)
    {
        InfoInFunction << "Could not find " << objName << endl;
    }
    return nullptr;
}


Foam::IOobjectList Foam::IOobjectList::lookup(const wordRe& matchName) const
{
    return lookupImpl(*this, matchName);
}


Foam::IOobjectList Foam::IOobjectList::lookup(const wordRes& matchName) const
{
    return lookupImpl(*this, matchName);
}


Foam::IOobjectList Foam::IOobjectList::lookup
(
    const wordHashSet& matchName
) const
{
    return lookupImpl(*this, matchName);
}


Foam::IOobjectList Foam::IOobjectList::lookupClass(const word& clsName) const
{
    return lookupClassImpl(*this, clsName, predicates::always());
}


Foam::HashTable<Foam::wordHashSet> Foam::IOobjectList::classes() const
{
    return classesImpl(*this, predicates::always());
}


Foam::HashTable<Foam::wordHashSet>
Foam::IOobjectList::classes(const wordRe& matchName) const
{
    return classesImpl(*this, matchName);
}


Foam::HashTable<Foam::wordHashSet>
Foam::IOobjectList::classes(const wordRes& matchName) const
{
    return classesImpl(*this, matchName);
}


Foam::HashTable<Foam::wordHashSet>
Foam::IOobjectList::classes(const wordHashSet& matchName) const
{
    return classesImpl(*this, matchName);
}


Foam::wordList Foam::IOobjectList::names() const
{
    return HashPtrTable<IOobject>::toc();
}


Foam::wordList Foam::IOobjectList::sortedNames() const
{
    return HashPtrTable<IOobject>::sortedToc();
}


Foam::wordList Foam::IOobjectList::names(const bool syncPar) const
{
    wordList objNames(HashPtrTable<IOobject>::toc());

    checkNames(objNames, syncPar);
    return objNames;
}


Foam::wordList Foam::IOobjectList::names
(
    const word& clsName
) const
{
    return namesImpl(*this, clsName, predicates::always(), false);
}


Foam::wordList Foam::IOobjectList::sortedNames
(
    const word& clsName
) const
{
    return namesImpl(*this, clsName, predicates::always(), true);
}


Foam::wordList Foam::IOobjectList::names
(
    const word& clsName,
    const bool syncPar
) const
{
    wordList objNames(namesImpl(*this, clsName, predicates::always(), false));

    checkNames(objNames, syncPar);
    return objNames;
}


Foam::wordList Foam::IOobjectList::names
(
    const word& clsName,
    const wordRe& matchName
) const
{
    return namesImpl(*this, clsName, matchName, false);
}


Foam::wordList Foam::IOobjectList::sortedNames
(
    const word& clsName,
    const wordRe& matchName
) const
{
    return namesImpl(*this, clsName, matchName, true);
}


Foam::wordList Foam::IOobjectList::names
(
    const word& clsName,
    const wordRe& matchName,
    const bool syncPar
) const
{
    wordList objNames(namesImpl(*this, clsName, matchName, false));

    checkNames(objNames, syncPar);
    return objNames;
}


Foam::wordList Foam::IOobjectList::names
(
    const word& clsName,
    const wordRes& matchName
) const
{
    return namesImpl(*this, clsName, matchName, false);
}


Foam::wordList Foam::IOobjectList::sortedNames
(
    const word& clsName,
    const wordRes& matchName
) const
{
    return namesImpl(*this, clsName, matchName, true);
}


Foam::wordList Foam::IOobjectList::names
(
    const word& clsName,
    const wordRes& matchName,
    const bool syncPar
) const
{
    wordList objNames(namesImpl(*this, clsName, matchName, false));

    checkNames(objNames, syncPar);
    return objNames;
}


Foam::wordList Foam::IOobjectList::names
(
    const word& clsName,
    const wordHashSet& matchName
) const
{
    return namesImpl(*this, clsName, matchName, false);
}


Foam::wordList Foam::IOobjectList::sortedNames
(
    const word& clsName,
    const wordHashSet& matchName
) const
{
    return namesImpl(*this, clsName, matchName, true);
}


Foam::wordList Foam::IOobjectList::names
(
    const word& clsName,
    const wordHashSet& matchName,
    const bool syncPar
) const
{
    wordList objNames(namesImpl(*this, clsName, matchName, false));

    checkNames(objNames, syncPar);
    return objNames;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::IOobjectList::operator=(IOobjectList&& list)
{
    transfer(list);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const IOobjectList& list)
{
    os << nl << list.size() << nl << token::BEGIN_LIST << nl;

    forAllConstIters(list, it)
    {
        os << it.key() << token::SPACE << it.object()->headerClassName() << nl;
    }

    os << token::END_LIST;
    os.check(FUNCTION_NAME);

    return os;
}


// ************************************************************************* //
