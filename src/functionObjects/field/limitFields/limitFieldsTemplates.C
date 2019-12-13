/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

#include "limitFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::limitFields::limitField(const word& fieldName)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    auto* fieldPtr = getObjectPtr<VolFieldType>(fieldName);
    if (!fieldPtr)
    {
        return false;
    }

    auto& field = *fieldPtr;

    Log << "    Limiting field " << fieldName << ":";

    const dimensionedScalar eps("eps", field.dimensions(), ROOTVSMALL);

    if (limit_ & MIN)
    {
        volScalarField mField(typeName + ":mag" + field.name(), mag(field));
        Log << " min(|" << gMin(mField) << "|)";
        //field.normalise();
        field /= mag(field) + eps;
        mField.max(dimensionedScalar("min", field.dimensions(), min_));
        field *= mField;
    }

    if (limit_ & MAX)
    {
        volScalarField mField(typeName + ":mag" + field.name(), mag(field));
        Log << " max(|" << gMax(mField) << "|)";
        //field.normalise();
        field /= mag(field) + eps;
        mField.min(dimensionedScalar("max", field.dimensions(), max_));
        field *= mField;
    }

    return true;
}


// ************************************************************************* //
