/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
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

Description
    Class for handling debugging switches.

\*---------------------------------------------------------------------------*/

#include "debug.H"
#include "dictionary.H"
#include "IFstream.H"
#include "OSspecific.H"
#include "Ostream.H"
//#include "HashTable.H"
#include "demandDrivenData.H"
#include "simpleObjectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace debug
{

//! \cond ignoreDocumentation - local scope
dictionary* controlDictPtr_(NULL);
dictionary* debugSwitchesPtr_(NULL);
dictionary* infoSwitchesPtr_(NULL);
dictionary* optimisationSwitchesPtr_(NULL);

// to ensure controlDictPtr_ is deleted at the end of the run
class deleteControlDictPtr
{
public:

    deleteControlDictPtr()
    {}

    ~deleteControlDictPtr()
    {
        deleteDemandDrivenData(controlDictPtr_);
    }
};

deleteControlDictPtr deleteControlDictPtr_;


// Debug switch read and write callback tables.
simpleObjectRegistry* debugObjectsPtr_(NULL);
simpleObjectRegistry* infoObjectsPtr_(NULL);
simpleObjectRegistry* optimisationObjectsPtr_(NULL);
class deleteDebugSwitchPtr
{
public:

    deleteDebugSwitchPtr()
    {}

    ~deleteDebugSwitchPtr()
    {
        deleteDemandDrivenData(debugObjectsPtr_);
        deleteDemandDrivenData(infoObjectsPtr_);
        deleteDemandDrivenData(optimisationObjectsPtr_);

    }
};

deleteDebugSwitchPtr deleteDebugSwitchPtr_;
//! \endcond


} // End namespace debug
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::dictionary& Foam::debug::controlDict()
{
    if (!controlDictPtr_)
    {
        fileNameList controlDictFiles = findEtcFiles("controlDict", true);
        controlDictPtr_ = new dictionary();
        forAllReverse(controlDictFiles, cdfi)
        {
            controlDictPtr_->merge
            (
                dictionary(IFstream(controlDictFiles[cdfi])())
            );
        }
    }

    return *controlDictPtr_;
}


Foam::dictionary& Foam::debug::switchSet
(
    const char* subDictName,
    dictionary*& subDictPtr
)
{
    if (!subDictPtr)
    {
        entry* ePtr = controlDict().lookupEntryPtr
        (
            subDictName, false, false
        );

        if (!ePtr || !ePtr->isDict())
        {
            cerr<< "debug::switchSet(const char*, dictionary*&):\n"
                << "    Cannot find " <<  subDictName << " in dictionary "
                << controlDict().name().c_str()
                << std::endl << std::endl;

            ::exit(1);
        }

        subDictPtr = &ePtr->dict();
    }

    return *subDictPtr;
}


Foam::dictionary& Foam::debug::debugSwitches()
{
    return switchSet("DebugSwitches", debugSwitchesPtr_);
}


Foam::dictionary& Foam::debug::infoSwitches()
{
    return switchSet("InfoSwitches", infoSwitchesPtr_);
}


Foam::dictionary& Foam::debug::optimisationSwitches()
{
    return switchSet("OptimisationSwitches", optimisationSwitchesPtr_);
}


int Foam::debug::debugSwitch(const char* name, const int defaultValue)
{
    return debugSwitches().lookupOrAddDefault
    (
        name, defaultValue, false, false
    );
}


int Foam::debug::infoSwitch(const char* name, const int defaultValue)
{
    return infoSwitches().lookupOrAddDefault
    (
        name, defaultValue, false, false
    );
}


int Foam::debug::optimisationSwitch(const char* name, const int defaultValue)
{
    return optimisationSwitches().lookupOrAddDefault
    (
        name, defaultValue, false, false
    );
}


void Foam::debug::addDebugObject(const char* name, simpleRegIOobject* obj)
{
    debugObjects().insert(name, obj);
}


void Foam::debug::addInfoObject(const char* name, simpleRegIOobject* obj)
{
    infoObjects().insert(name, obj);
}


void Foam::debug::addOptimisationObject
(
    const char* name,
    simpleRegIOobject* obj
)
{
    optimisationObjects().insert(name, obj);
}


Foam::simpleObjectRegistry& Foam::debug::debugObjects()
{
    if (!debugObjectsPtr_)
    {
        debugObjectsPtr_ = new simpleObjectRegistry(1000);
    }

    return *debugObjectsPtr_;
}


Foam::simpleObjectRegistry& Foam::debug::infoObjects()
{
    if (!infoObjectsPtr_)
    {
        infoObjectsPtr_ = new simpleObjectRegistry(1000);
    }

    return *infoObjectsPtr_;
}


Foam::simpleObjectRegistry& Foam::debug::optimisationObjects()
{
    if (!optimisationObjectsPtr_)
    {
        optimisationObjectsPtr_ = new simpleObjectRegistry(1000);
    }

    return *optimisationObjectsPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
