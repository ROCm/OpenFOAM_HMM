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

#include "mappedFixedInternalValueFvPatchField.H"
#include "UIndirectList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::mappedFixedInternalValueFvPatchField<Type>::
mappedFixedInternalValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mappedFixedValueFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::mappedFixedInternalValueFvPatchField<Type>::
mappedFixedInternalValueFvPatchField
(
    const mappedFixedInternalValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mappedFixedValueFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::mappedFixedInternalValueFvPatchField<Type>::
mappedFixedInternalValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mappedFixedValueFvPatchField<Type>(p, iF, dict)
{}


template<class Type>
Foam::mappedFixedInternalValueFvPatchField<Type>::
mappedFixedInternalValueFvPatchField
(
    const mappedFixedInternalValueFvPatchField<Type>& ptf
)
:
    mappedFixedValueFvPatchField<Type>(ptf)
{}


template<class Type>
Foam::mappedFixedInternalValueFvPatchField<Type>::
mappedFixedInternalValueFvPatchField
(
    const mappedFixedInternalValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mappedFixedValueFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::mappedFixedInternalValueFvPatchField<Type>::updateCoeffs()
{
    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;

    if (this->updated())
    {
        return;
    }

    // Retrieve the neighbour values and assign to this patch boundary field
    mappedFixedValueFvPatchField<Type>::updateCoeffs();

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(this->patch().patch());
    const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
    const label samplePatchI = mpp.samplePolyPatch().index();
    const fvPatch& nbrPatch = nbrMesh.boundary()[samplePatchI];


    // Retrieve the neighbour field
    const fvPatchField<Type>& nbrField =
        nbrPatch.template lookupPatchField<FieldType, Type>
        (
            this->dimensionedInternalField().name()
        );

    // Retrieve the neighbour patch internal field
    Field<Type> nbrIntFld(nbrField.patchInternalField());

    switch (mpp.mode())
    {
        case (mappedPatchBase::NEARESTPATCHFACEAMI):
        {
            // Retrieve the neighbour patch internal field
            nbrIntFld = mpp.AMI().interpolateToSource(nbrIntFld);
            break;
        }
        default:
        {
            const mapDistribute& distMap = mpp.map();

            mapDistribute::distribute
            (
                Pstream::defaultCommsType,
                distMap.schedule(),
                distMap.constructSize(),
                distMap.subMap(),           // what to send
                distMap.constructMap(),     // what to receive
                nbrIntFld
            );
            break;
        }
    }

    // Assign (this) patch internal field to its neighbour values
    Field<Type>& intFld = const_cast<Field<Type>&>(this->internalField());
    UIndirectList<Type>(intFld, this->patch().faceCells()) = nbrIntFld;
}


template<class Type>
void Foam::mappedFixedInternalValueFvPatchField<Type>::write
(
    Ostream& os
) const
{
    mappedFixedValueFvPatchField<Type>::write(os);
}


// ************************************************************************* //
