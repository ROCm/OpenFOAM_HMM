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

\*---------------------------------------------------------------------------*/

#include "functionObjectList.H"
#include "Time.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::functionObject* Foam::functionObjectList::remove(const word& key)
{
    functionObject* ptr = 0;

    // Find index of existing functionObject
    HashTable<label>::iterator fnd = indices_.find(key);

    if (fnd != indices_.end())
    {
        // remove the pointer from the old list
        ptr = functions_.set(fnd(), 0).ptr();
        indices_.erase(fnd);
    }

    return ptr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjectList::functionObjectList
(
    const Time& t,
    const bool execution
)
:
    functions_(),
    indices_(),
    time_(t),
    parentDict_(t.controlDict()),
    execution_(execution),
    updated_(false)
{}


Foam::functionObjectList::functionObjectList
(
    const Time& t,
    const dictionary& parentDict,
    const bool execution
)
:
    functions_(),
    indices_(),
    time_(t),
    parentDict_(parentDict),
    execution_(execution),
    updated_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjectList::~functionObjectList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjectList::start()
{
    return read();
}


bool Foam::functionObjectList::execute()
{
    bool ok = true;

    if (execution_)
    {
        if (!updated_)
        {
            read();
        }

        forAllIter(PtrList<functionObject>, functions_, iter)
        {
            ok = iter().execute() && ok;
        }
    }

    return ok;
}


void Foam::functionObjectList::on()
{
    execution_ = true;
}


void Foam::functionObjectList::off()
{
    execution_ = false;
}


bool Foam::functionObjectList::read()
{
    bool ok = true;
    updated_ = execution_;

    // avoid reading/initializing if execution is off
    if (!execution_)
    {
        return ok;
    }

    // Update existing and add new functionObjects
    const entry* entryPtr = parentDict_.lookupEntryPtr("functions",false,false);
    if (entryPtr)
    {
        PtrList<functionObject> newPtrs;
        HashTable<label> newIndices;

        label nFunc = 0;

        if (entryPtr->isDict())
        {
            // a dictionary of functionObjects
            const dictionary& functionDicts = entryPtr->dict();
            newPtrs.setSize(functionDicts.size());

            forAllConstIter(dictionary, functionDicts, iter)
            {
                // safety:
                if (!iter().isDict())
                {
                    continue;
                }
                const word& key = iter().keyword();
                const dictionary& dict = iter().dict();

                functionObject* objPtr = remove(key);
                if (objPtr)
                {
                    // existing functionObject
                    ok = objPtr->read(dict) && ok;
                }
                else
                {
                    // new functionObject
                    objPtr = functionObject::New(key, time_, dict).ptr();
                    ok = objPtr->start() && ok;
                }

                newPtrs.set(nFunc, objPtr);
                newIndices.insert(key, nFunc);
                nFunc++;
            }
        }
        else
        {
            // a list of functionObjects
            PtrList<entry> functionDicts(entryPtr->stream());
            newPtrs.setSize(functionDicts.size());

            forAllIter(PtrList<entry>, functionDicts, iter)
            {
                // safety:
                if (!iter().isDict())
                {
                    continue;
                }
                const word& key = iter().keyword();
                const dictionary& dict = iter().dict();

                functionObject* objPtr = remove(key);
                if (objPtr)
                {
                    // existing functionObject
                    ok = objPtr->read(dict) && ok;
                }
                else
                {
                    // new functionObject
                    objPtr = functionObject::New(key, time_, dict).ptr();
                    ok = objPtr->start() && ok;
                }

                newPtrs.set(nFunc, objPtr);
                newIndices.insert(key, nFunc);
                nFunc++;
            }
        }

        // safety:
        newPtrs.setSize(nFunc);

        // update PtrList of functionObjects
        // also deletes existing, unused functionObjects
        functions_.transfer(newPtrs);
        indices_.transfer(newIndices);
    }
    else
    {
        functions_.clear();
        indices_.clear();
    }

    return ok;
}


// ************************************************************************* //
