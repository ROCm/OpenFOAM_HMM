/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

template <class Type>
void Foam::explicitSource::addSource
(
    fvMatrix<Type>& Eqn,
    const Type& sourceData
) const
{
    Type data = sourceData;
    if (volumeMode_ == vmAbsolute)
    {
        // Convert to specific quantity
        data /= V_;
    }

    DimensionedField<Type, volMesh> rhs
    (
        IOobject
        (
            "rhs",
            Eqn.psi().mesh().time().timeName(),
            Eqn.psi().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        Eqn.psi().mesh(),
        dimensioned<Type>
        (
            "zero",
            Eqn.dimensions()/dimVolume,
            pTraits<Type>::zero
        )
    );
    UIndirectList<Type>(rhs, this->cells()) = data;    

    Eqn -= rhs;
}


template <class Type>
void Foam::explicitSource::addField
(
    HashTable<Type>& fields,
    const wordList& fieldTypes,
    const wordList& fieldNames,
    const dictionary& fieldDataDict
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> geometricField;

    forAll (fieldTypes, fieldI)
    {
        word fieldName = fieldNames[fieldI];
        word fieldType = fieldTypes[fieldI];

        if
        (
            (
                fieldType
             == GeometricField<Type, fvPatchField, volMesh>::typeName
            ) &&
            (
                this->mesh().foundObject<geometricField>(fieldName)
            )
        )
        {
            Type fieldValue = fieldDataDict.lookupOrDefault<Type>
            (
                fieldName,
                pTraits<Type>::zero
            );

            fields.insert(fieldName, fieldValue);
        }
    }
}


// ************************************************************************* //
