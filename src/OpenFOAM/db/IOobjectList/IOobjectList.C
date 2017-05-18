/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "IOobjectList.H"
#include "Time.H"
#include "OSspecific.H"
#include "IOList.H"
#include "stringListOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IOobjectList::IOobjectList(const label nIoObjects)
:
    HashPtrTable<IOobject>(nIoObjects)
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
    word newInstance = instance;

    if (!isDir(db.path(instance)))
    {
        newInstance = db.time().findInstancePath(instant(instance));

        if (newInstance.empty())
        {
            return;
        }
    }

    // Create a list of file names in this directory
    fileNameList ObjectNames =
        readDir(db.path(newInstance, db.dbDir()/local), fileName::FILE);

    forAll(ObjectNames, i)
    {
        IOobject* objectPtr = new IOobject
        (
            ObjectNames[i],
            newInstance,
            local,
            db,
            r,
            w,
            registerObject
        );

        // Use object with local scope
        if (objectPtr->typeHeaderOk<IOList<label>>(false))
        {
            insert(ObjectNames[i], objectPtr);
        }
        else
        {
            delete objectPtr;
        }
    }
}


Foam::IOobjectList::IOobjectList(const IOobjectList& iolist)
:
    HashPtrTable<IOobject>(iolist)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::IOobjectList::~IOobjectList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::IOobjectList::add(IOobject& io)
{
    return insert(io.name(), &io);
}


bool Foam::IOobjectList::remove(IOobject& io)
{
    return erase(io.name());
}


Foam::IOobject* Foam::IOobjectList::lookup(const word& name) const
{
    HashPtrTable<IOobject>::const_iterator iter = find(name);

    if (iter.found())
    {
        if (IOobject::debug)
        {
            InfoInFunction << "Found " << name << endl;
        }

        return const_cast<IOobject*>(*iter);
    }
    else
    {
        if (IOobject::debug)
        {
            InfoInFunction << "Could not find " << name << endl;
        }

        return nullptr;
    }
}


Foam::IOobjectList Foam::IOobjectList::lookup(const wordRe& matcher) const
{
    IOobjectList results(size());

    forAllConstIters(*this, iter)
    {
        if (matcher.match(iter.key()))
        {
            if (IOobject::debug)
            {
                InfoInFunction << "Found " << iter.key() << endl;
            }

            results.insert
            (
                iter.key(),
                new IOobject(*(iter.object()))
            );
        }
    }

    return results;
}


Foam::IOobjectList Foam::IOobjectList::lookup(const wordReList& matcher) const
{
    wordReListMatcher mat(matcher);

    IOobjectList results(size());

    forAllConstIters(*this, iter)
    {
        if (mat.match(iter.key()))
        {
            if (IOobject::debug)
            {
                InfoInFunction << "Found " << iter.key() << endl;
            }

            results.insert
            (
                iter.key(),
                new IOobject(*(iter.object()))
            );
        }
    }

    return results;
}


Foam::IOobjectList Foam::IOobjectList::lookupClass(const word& clsName) const
{
    IOobjectList results(size());

    forAllConstIters(*this, iter)
    {
        if (iter()->headerClassName() == clsName)
        {
            if (IOobject::debug)
            {
                InfoInFunction << "Found " << iter.key() << endl;
            }

            results.insert
            (
                iter.key(),
                new IOobject(*(iter.object()))
            );
        }
    }

    return results;
}


Foam::wordList Foam::IOobjectList::names() const
{
    return HashPtrTable<IOobject>::toc();
}


Foam::wordList Foam::IOobjectList::sortedNames() const
{
    return HashPtrTable<IOobject>::sortedToc();
}


Foam::wordList Foam::IOobjectList::names
(
    const word& clsName
) const
{
    wordList objNames(size());

    label count = 0;
    forAllConstIters(*this, iter)
    {
        if (iter()->headerClassName() == clsName)
        {
            objNames[count++] = iter.key();
        }
    }

    objNames.setSize(count);

    return objNames;
}


Foam::wordList Foam::IOobjectList::names
(
    const word& clsName,
    const wordRe& matcher
) const
{
    wordList objNames = names(clsName);

    return wordList(objNames, findStrings(matcher, objNames));
}


Foam::wordList Foam::IOobjectList::names
(
    const word& clsName,
    const wordReList& matcher
) const
{
    wordList objNames = names(clsName);

    return wordList(objNames, findStrings(matcher, objNames));
}


Foam::wordList Foam::IOobjectList::sortedNames
(
    const word& clsName
) const
{
    wordList sortedLst = names(clsName);
    sort(sortedLst);

    return sortedLst;
}


Foam::wordList Foam::IOobjectList::sortedNames
(
    const word& clsName,
    const wordRe& matcher
) const
{
    wordList sortedLst = names(clsName, matcher);
    sort(sortedLst);

    return sortedLst;
}


Foam::wordList Foam::IOobjectList::sortedNames
(
    const word& clsName,
    const wordReList& matcher
) const
{
    wordList sortedLst = names(clsName, matcher);
    sort(sortedLst);

    return sortedLst;
}


// ************************************************************************* //
