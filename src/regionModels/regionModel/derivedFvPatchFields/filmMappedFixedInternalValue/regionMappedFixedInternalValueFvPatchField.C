/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "regionMappedFixedInternalValueFvPatchField.H"
#include "UIndirectList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::regionMappedFixedInternalValueFvPatchField<Type>::
regionMappedFixedInternalValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::regionMappedFixedInternalValueFvPatchField<Type>::
regionMappedFixedInternalValueFvPatchField
(
    const regionMappedFixedInternalValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::regionMappedFixedInternalValueFvPatchField<Type>::
regionMappedFixedInternalValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict)
{}


template<class Type>
Foam::regionMappedFixedInternalValueFvPatchField<Type>::
regionMappedFixedInternalValueFvPatchField
(
    const regionMappedFixedInternalValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf)
{}


template<class Type>
Foam::regionMappedFixedInternalValueFvPatchField<Type>::
regionMappedFixedInternalValueFvPatchField
(
    const regionMappedFixedInternalValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::regionMappedFixedInternalValueFvPatchField<Type>::updateCoeffs()
{
    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;

    if (this->updated())
    {
        return;
    }

    const regionModels::regionModel& region =
        this->db.lookupObject<regionModels::regionModel>
        (
            "surfaceFilmProperties"
        );

    const label regionPatchI = region.regionPatchID(this->patch().index())
    const directMappedPatchBase& mpp = region.mappedPatches()[patchI];

    // Force recalculation of mapping and schedule
    const mapDistribute& distMap = mpp.map();

    const polyMesh& nbrMesh = mpp.sampleMesh();
    const label samplePatchI = mpp.samplePolyPatch().index();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[samplePatchI];


    // Retrieve the neighbour field
    Field<Type> nbrField =
        nbrPatch.lookupPatchField<FieldType, Type>(fieldName_);

    mapDistribute::distribute
    (
        Pstream::defaultCommsType,
        distMap.schedule(),
        distMap.constructSize(),
        distMap.subMap(),           // what to send
        distMap.constructMap(),     // what to receive
        nbrField
    );

    this->operator==(nbrField);

    // Retrieve the neighbour patch internal field
    Field<Type> nbrIntField = nbrField.patchInternalField();
    mapDistribute::distribute
    (
        Pstream::defaultCommsType,
        distMap.schedule(),
        distMap.constructSize(),
        distMap.subMap(),           // what to send
        distMap.constructMap(),     // what to receive
        nbrIntField
    );

    // Assign (this) patch internal field to its neighbour values
    Field<Type>& intField = const_cast<Field<Type>&>(this->internalField());
    UIndirectList<Type>(intFld, this->patch().faceCells()) = nbrIntField;

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::regionMappedFixedInternalValueFvPatchField<Type>::write
(
    Ostream& os
) const
{
    fixedValueFvPatchField<Type>::write(os);
}


// ************************************************************************* //
