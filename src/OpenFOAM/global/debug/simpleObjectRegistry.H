/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Class
    Foam::simpleObjectRegistry

Description
    Object registry for simpleRegIOobject. Maintains ordering.

SourceFiles
    simpleObjectRegistry.C

\*---------------------------------------------------------------------------*/

#ifndef simpleObjectRegistry_H
#define simpleObjectRegistry_H

#include "Dictionary.H"
#include "simpleRegIOobject.H"
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class dictionary;

/*---------------------------------------------------------------------------*\
                   Class simpleObjectRegistryEntry Declaration
\*---------------------------------------------------------------------------*/

class simpleObjectRegistryEntry
:
    public Dictionary<simpleObjectRegistryEntry>::link,
    public List<simpleRegIOobject*>
{
public:

    //- Construct with a single object (list size == 1)
    explicit simpleObjectRegistryEntry(simpleRegIOobject* obj)
    :
        List<simpleRegIOobject*>(label(1), obj)
    {}

    //- Construct with a List of objects
    explicit simpleObjectRegistryEntry(const List<simpleRegIOobject*>& data)
    :
        List<simpleRegIOobject*>(data)
    {}
};


/*---------------------------------------------------------------------------*\
                      Class simpleObjectRegistry Declaration
\*---------------------------------------------------------------------------*/

class simpleObjectRegistry
:
    public Dictionary<simpleObjectRegistryEntry>
{
public:

    // Constructors

        //- Construct with given or default (128) table capacity
        explicit simpleObjectRegistry(const label size = 128)
        :
            Dictionary<simpleObjectRegistryEntry>(size)
        {}


    // Member Functions

        //- Set values (invoke callbacks) from dictionary entries
        //  Reporting honours the infoDetailLevel
        void setValues(const dictionary& dict, bool report=false);

        //- Set named value, but also handle embedded name=value syntax
        //  Reporting honours the infoDetailLevel
        void setNamedInt(std::string name, int val, bool report=false);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
