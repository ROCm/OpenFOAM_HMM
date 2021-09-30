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
#include "ITstream.H"
#include "StringStream.H"
#include "int.H"
#include "floatScalar.H"

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


void Foam::simpleObjectRegistry::setNamedValue
(
    std::string name,
    int val,
    bool report
)
{
    // Report enables output, but respect DetailInfo state as well.
    // The local log variable captures this logic.

    const bool log = (report && Foam::infoDetailLevel > 0);

    token tok(static_cast<label>(val));

    // Handle name=value
    const auto eq = name.find('=');

    if (eq != std::string::npos)
    {
        std::string strval(name.substr(eq+1));
        name.erase(eq);  // Truncate the name

        float fvalue(val);

        if (Foam::readInt(strval, val))
        {
            tok = static_cast<label>(val);
        }
        else if (Foam::readFloat(strval, fvalue))
        {
            tok = fvalue;
        }
        // Treat 'name=' like 'name' (ie, no value parameter)
        // silently ignore 'name=junk', but could warn
    }


    simpleObjectRegistryEntry* objPtr = this->find(name.c_str());

    if (objPtr)
    {
        // The generic interface requires an Istream.
        ITstream is("", tokenList(Foam::one{}, tok));

        Log << name.c_str() << '=' << tok << nl;

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
