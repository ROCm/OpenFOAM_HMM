/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "uniformFixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformFixedValueFvPatchField<Type>::uniformFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    uniformValue_(nullptr)
{}


template<class Type>
Foam::uniformFixedValueFvPatchField<Type>::uniformFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& fld
)
:
    fixedValueFvPatchField<Type>(p, iF, fld),
    uniformValue_(nullptr)
{}


template<class Type>
Foam::uniformFixedValueFvPatchField<Type>::uniformFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict, false),
    uniformValue_(PatchFunction1<Type>::New(p.patch(), "uniformValue", dict))
{
    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        this->evaluate();
    }
}


template<class Type>
Foam::uniformFixedValueFvPatchField<Type>::uniformFixedValueFvPatchField
(
    const uniformFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(p, iF),   // Don't map
    uniformValue_(ptf.uniformValue_.clone(p.patch()))
{
    if (mapper.direct() && !mapper.hasUnmapped())
    {
        // Use mapping instead of re-evaluation
        this->map(ptf, mapper);
    }
    else
    {
        // Evaluate since value not mapped
        this->evaluate();
    }
}


template<class Type>
Foam::uniformFixedValueFvPatchField<Type>::uniformFixedValueFvPatchField
(
    const uniformFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    uniformValue_(ptf.uniformValue_.clone(this->patch().patch()))
{}


template<class Type>
Foam::uniformFixedValueFvPatchField<Type>::uniformFixedValueFvPatchField
(
    const uniformFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    uniformValue_(ptf.uniformValue_.clone(this->patch().patch()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformFixedValueFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchField<Type>::autoMap(mapper);
    uniformValue_().autoMap(mapper);

    if (uniformValue_().constant())
    {
        // If mapper is not dependent on time we're ok to evaluate
        this->evaluate();
    }
}


template<class Type>
void Foam::uniformFixedValueFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const uniformFixedValueFvPatchField& tiptf =
        refCast<const uniformFixedValueFvPatchField>(ptf);

    uniformValue_().rmap(tiptf.uniformValue_(), addr);
}


template<class Type>
void Foam::uniformFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();
    fvPatchField<Type>::operator==(uniformValue_->value(t));
    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::uniformFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    uniformValue_->writeData(os);
    this->writeEntry("value", os);
}


// ************************************************************************* //
