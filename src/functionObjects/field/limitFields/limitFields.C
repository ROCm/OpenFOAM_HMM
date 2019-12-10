/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd
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

#include "limitFields.H"
#include "fieldTypes.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(limitFields, 0);
    addToRunTimeSelectionTable(functionObject, limitFields, dictionary);
}
}

const Foam::Enum
<
    Foam::functionObjects::limitFields::limitType
>
Foam::functionObjects::limitFields::limitTypeNames_
({
    { limitType::MIN, "min" },
    { limitType::MAX, "max" },
    { limitType::BOTH, "both" },
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::limitFields::limitScalarField
(
    const word& fieldName
)
{
    auto* fieldPtr = obr_.getObjectPtr<volScalarField>(fieldName);
    if (!fieldPtr)
    {
        return false;
    }

    auto& field = *fieldPtr;

    if (limit_ & MIN)
    {
        Log << ": min(" << gMin(field) << ")";
        field.max(dimensionedScalar("", field.dimensions(), min_));
    }

    if (limit_ & MAX)
    {
        Log << ": max(" << gMax(field) << ")";
        field.min(dimensionedScalar("", field.dimensions(), max_));
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::limitFields::limitFields
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    limit_(MIN),
    fieldSet_(mesh_),
    min_(-VGREAT),
    max_(VGREAT)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::limitFields::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {
        Info<< type() << " " << name() << ":" << nl;

        limit_ = limitTypeNames_.get("limit", dict);

        if (limit_ & MIN)
        {
            min_ = dict.get<scalar>("min");
            Info<< "    Imposing lower limit " << min_ << nl;
        }

        if (limit_ & MAX)
        {
            max_ = dict.get<scalar>("max");
            Info<< "    Imposing upper limit " << max_ << nl;
        }

        fieldSet_.read(dict);

        Info<< endl;

        return true;
    }

    return false;
}


bool Foam::functionObjects::limitFields::execute()
{
    fieldSet_.updateSelection();

    Log << type() << " " << name() << ":" << nl;

    label count = 0;
    for (const word& fieldName : fieldSet_.selectionNames())
    {
        if
        (
            limitScalarField(fieldName)
         || limitField<vector>(fieldName)
         || limitField<sphericalTensor>(fieldName)
         || limitField<symmTensor>(fieldName)
         || limitField<tensor>(fieldName)
        )
        {
            ++count;
        }
    }

    if (debug)
    {
        Log << " - limited " << count << '/'
            << fieldSet_.selectionNames().size() << " fields";
    }

    Log << endl;

    return true;
}


bool Foam::functionObjects::limitFields::write()
{
    for (const word& fieldName : fieldSet_.selectionNames())
    {
        lookupObject<regIOobject>(fieldName).write();
    }

    return true;
}


// ************************************************************************* //
