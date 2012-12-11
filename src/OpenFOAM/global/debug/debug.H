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

Namespace
    Foam::debug

Description
    Namespace for handling debugging switches.

SourceFiles
    debug.C

\*---------------------------------------------------------------------------*/

#ifndef debug_H
#define debug_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class dictionary;
class Istream;
class Ostream;
class simpleRegIOobject;
class simpleObjectRegistry;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace debug
{
    //- The central control dictionary.
    //  Located in ~/.OpenFOAM/VERSION or $WM_PROJECT_DIR/etc
    //  \sa Foam::findEtcFile()
    dictionary& controlDict();

    //- The DebugSwitches sub-dictionary in the central controlDict.
    dictionary& debugSwitches();

    //- The InfoSwitches sub-dictionary in the central controlDict.
    dictionary& infoSwitches();

    //- The OptimisationSwitches sub-dictionary in the central controlDict.
    dictionary& optimisationSwitches();

    //- Lookup debug switch or add default value.
    int debugSwitch(const char* name, const int defaultValue=0);

    //- Lookup info switch or add default value.
    int infoSwitch(const char* name, const int defaultValue=0);

    //- Lookup optimisation switch or add default value.
    int optimisationSwitch(const char* name, const int defaultValue=0);

    //- Internal function to lookup a sub-dictionary from controlDict.
    dictionary& switchSet(const char* subDictName, dictionary*& subDictPtr);



    // Registered debug switches

        //- Register debug switch read/write object
        void addDebugObject(const char* name, simpleRegIOobject* obj);

        //- Register debug switch read/write object
        void addInfoObject(const char* name, simpleRegIOobject* obj);

        //- Register debug switch read/write object
        void addOptimisationObject(const char* name, simpleRegIOobject* obj);


        //- Get access to registered debug switch objects
        simpleObjectRegistry& debugObjects();

        //- Get access to registered debug switch objects
        simpleObjectRegistry& infoObjects();

        //- Get access to registered debug switch objects
        simpleObjectRegistry& optimisationObjects();


} // End namespace debug


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
