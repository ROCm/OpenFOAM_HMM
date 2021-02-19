/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "simpleObjectRegistry.H"
#include "dictionary.H"
#include "StringStream.H"
#include "int.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::simpleObjectRegistry::setValues
(
    const dictionary& dict,
    bool report
)
{
    // Report enables output, but respect DetailInfo state as well.
    // The local log variable captures this logic.

    const bool log = (report && Foam::infoDetailLevel > 0);

    for (const entry& dEntry : dict)
    {
        const word& name = dEntry.keyword();

        simpleObjectRegistryEntry* objPtr = this->find(name);

        if (objPtr)
        {
            Log << "    " << dEntry << nl;

            const List<simpleRegIOobject*>& objects = *objPtr;

            if (dEntry.isDict())
            {
                OStringStream os;
                os  << dEntry.dict();
                IStringStream is(os.str());

                // Or alternatively?
                // ITstream is(name, dEntry.dict().tokens());

                for (simpleRegIOobject* obj : objects)
                {
                    is.rewind();
                    obj->readData(is);
                }
            }
            else
            {
                for (simpleRegIOobject* obj : objects)
                {
                    obj->readData(dEntry.stream());
                }
            }
        }
        else
        {
            Log << "    " << name << " (unregistered)" << nl;
        }
    }
}


void Foam::simpleObjectRegistry::setNamedInt
(
    std::string name,
    int val,
    bool report
)
{
    // Report enables output, but respect DetailInfo state as well.
    // The local log variable captures this logic.

    const bool log = (report && Foam::infoDetailLevel > 0);


    // Handle name=value
    const auto eq = name.find('=');

    if (eq != std::string::npos)
    {
        int intval = 0;

        if (readInt(name.substr(eq+1), intval))
        {
            val = intval;
        }
        // Could warn about bad entry

        name.resize(eq);  // Truncate the name
    }


    simpleObjectRegistryEntry* objPtr = this->find(name.c_str());

    if (objPtr)
    {
        // The generic interface requires an Istream.
        IStringStream is(std::to_string(val));

        // Or alternatively?
        // ITstream is("input", tokenList(1, token(label(val))));

        Log << name.c_str() << '=' << val << nl;

        const List<simpleRegIOobject*>& objects = *objPtr;

        for (simpleRegIOobject* obj : objects)
        {
            is.rewind();
            obj->readData(is);
        }
    }
    else
    {
        Log << name.c_str() << " (unregistered)" << nl;
    }
}


// ************************************************************************* //
