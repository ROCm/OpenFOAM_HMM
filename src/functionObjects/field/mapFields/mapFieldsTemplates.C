/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "fvMesh.H"
#include "polyPatch.H"
#include "lduSchedule.H"
#include "meshToMesh.H"

template<class Type>
void Foam::functionObjects::mapFields::evaluateConstraintTypes
(
    GeometricField<Type, fvPatchField, volMesh>& fld
) const
{
    auto& fldBf = fld.boundaryFieldRef();

    const UPstream::commsTypes commsType(UPstream::defaultCommsType);

    if
    (
        commsType == UPstream::commsTypes::blocking
     || commsType == UPstream::commsTypes::nonBlocking
    )
    {
        const label startOfRequests = UPstream::nRequests();

        forAll(fldBf, patchi)
        {
            fvPatchField<Type>& tgtField = fldBf[patchi];

            if
            (
                tgtField.type() == tgtField.patch().patch().type()
             && polyPatch::constraintType(tgtField.patch().patch().type())
            )
            {
                tgtField.initEvaluate(commsType);
            }
        }

        // Wait for outstanding requests
        if
        (
            UPstream::parRun()
         && commsType == UPstream::commsTypes::nonBlocking
        )
        {
            UPstream::waitRequests(startOfRequests);
        }

        forAll(fldBf, patchi)
        {
            fvPatchField<Type>& tgtField = fldBf[patchi];

            if
            (
                tgtField.type() == tgtField.patch().patch().type()
             && polyPatch::constraintType(tgtField.patch().patch().type())
            )
            {
                tgtField.evaluate(commsType);
            }
        }
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule =
            fld.mesh().globalData().patchSchedule();

        for (const auto& schedEval : patchSchedule)
        {
            const label patchi = schedEval.patch;

            fvPatchField<Type>& tgtField = fldBf[patchi];

            if
            (
                tgtField.type() == tgtField.patch().patch().type()
             && polyPatch::constraintType(tgtField.patch().patch().type())
            )
            {
                if (schedEval.init)
                {
                    tgtField.initEvaluate(commsType);
                }
                else
                {
                    tgtField.evaluate(commsType);
                }
            }
        }
    }
}


template<class Type>
bool Foam::functionObjects::mapFields::mapFieldType() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    const fvMesh& mapRegion = mapRegionPtr_();

    wordList fieldNames(this->mesh_.names(VolFieldType::typeName));

    const labelList selected(fieldNames_.matching(fieldNames));

    for (const label fieldi : selected)
    {
        const word& fieldName = fieldNames[fieldi];

        const VolFieldType& field = lookupObject<VolFieldType>(fieldName);

        if (!mapRegion.foundObject<VolFieldType>(fieldName))
        {
            VolFieldType* tmappedField =
                new VolFieldType
                (
                    IOobject
                    (
                        fieldName,
                        time_.timeName(),
                        mapRegion,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mapRegion,
                    dimensioned<Type>(field.dimensions(), Zero)
                );

            tmappedField->store();
        }

        VolFieldType& mappedField =
            mapRegion.template lookupObjectRef<VolFieldType>(fieldName);

        mappedField = interpPtr_->mapTgtToSrc(field);

        Log << "    " << fieldName << ": interpolated";

        evaluateConstraintTypes(mappedField);
    }

    return !selected.empty();
}


template<class Type>
bool Foam::functionObjects::mapFields::writeFieldType() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    const fvMesh& mapRegion = mapRegionPtr_();

    wordList fieldNames(this->mesh_.names(VolFieldType::typeName));

    const labelList selected(fieldNames_.matching(fieldNames));

    for (const label fieldi : selected)
    {
        const word& fieldName = fieldNames[fieldi];

        const VolFieldType& mappedField =
            mapRegion.template lookupObject<VolFieldType>(fieldName);

        mappedField.write();

        Log << "    " << fieldName << ": written";
    }

    return !selected.empty();
}


// ************************************************************************* //
