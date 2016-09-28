/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include "readFields.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(readFields, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::readFields::readFields
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    fieldSet_(),
    log_(true)
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (isA<fvMesh>(obr_))
    {
        read(dict);

        // Fields should all be present from start time so read on construction
        execute();
    }
    else
    {
        active_ = false;
        WarningInFunction
            << "No fvMesh available, deactivating " << name_
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::readFields::~readFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::readFields::read(const dictionary& dict)
{
    if (active_)
    {
        log_.readIfPresent("log", dict);
        dict.lookup("fields") >> fieldSet_;
    }
}


void Foam::readFields::execute()
{
    if (active_)
    {
        if (log_) Info << type() << " " << name_ << ":" << nl;

        bool loaded = false;
        forAll(fieldSet_, fieldI)
        {
            const word& fieldName = fieldSet_[fieldI];

            // Load field if necessary
            loaded = loadField<scalar>(fieldName) || loaded;
            loaded = loadField<vector>(fieldName) || loaded;
            loaded = loadField<sphericalTensor>(fieldName) || loaded;
            loaded = loadField<symmTensor>(fieldName) || loaded;
            loaded = loadField<tensor>(fieldName) || loaded;
        }

        if (log_)
        {
            if (!loaded)
            {
                Info<< "    no fields loaded" << endl;
            }
            Info<< endl;
        }
    }
}


void Foam::readFields::end()
{
    // Do nothing
}


void Foam::readFields::timeSet()
{
    // Do nothing
}


void Foam::readFields::write()
{
    // Do nothing
}


// ************************************************************************* //
