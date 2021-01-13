/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "partialSlipFvPatchField.H"
#include "symmTransformField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::partialSlipFvPatchField<Type>::partialSlipFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    parent_bctype(p, iF),
    refValue_(p.size(), Zero),
    valueFraction_(p.size(), 1.0),
    writeValue_(false)
{}


template<class Type>
Foam::partialSlipFvPatchField<Type>::partialSlipFvPatchField
(
    const partialSlipFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(ptf, p, iF, mapper),
    refValue_(ptf.refValue_, mapper),
    valueFraction_(ptf.valueFraction_, mapper),
    writeValue_(ptf.writeValue_)
{}


template<class Type>
Foam::partialSlipFvPatchField<Type>::partialSlipFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF),
    refValue_(p.size(), Zero),
    valueFraction_("valueFraction", dict, p.size()),
    writeValue_(dict.getOrDefault("writeValue", false))
{
    this->patchType() = dict.getOrDefault<word>("patchType", word::null);

    // Backwards compatibility - leave refValue as zero unless specified
    if (dict.found("refValue"))
    {
        refValue_ = Field<Type>("refValue", dict, p.size());
    }

    evaluate();
}


template<class Type>
Foam::partialSlipFvPatchField<Type>::partialSlipFvPatchField
(
    const partialSlipFvPatchField<Type>& ptf
)
:
    parent_bctype(ptf),
    refValue_(ptf.refValue_),
    valueFraction_(ptf.valueFraction_),
    writeValue_(ptf.writeValue_)
{}


template<class Type>
Foam::partialSlipFvPatchField<Type>::partialSlipFvPatchField
(
    const partialSlipFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    parent_bctype(ptf, iF),
    refValue_(ptf.refValue_),
    valueFraction_(ptf.valueFraction_),
    writeValue_(ptf.writeValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::partialSlipFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    parent_bctype::autoMap(m);
    refValue_.autoMap(m);
    valueFraction_.autoMap(m);
}


template<class Type>
void Foam::partialSlipFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    parent_bctype::rmap(ptf, addr);

    const auto& dmptf =
        refCast<const partialSlipFvPatchField<Type>>(ptf);

    refValue_.rmap(dmptf.refValue_, addr);
    valueFraction_.rmap(dmptf.valueFraction_, addr);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::partialSlipFvPatchField<Type>::snGrad() const
{
    tmp<vectorField> nHat = this->patch().nf();
    const Field<Type> pif(this->patchInternalField());

    return
    (
        valueFraction_*refValue_
      + (1.0 - valueFraction_)*transform(I - sqr(nHat), pif) - pif
    )*this->patch().deltaCoeffs();
}


template<class Type>
void Foam::partialSlipFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    tmp<vectorField> nHat = this->patch().nf();

    Field<Type>::operator=
    (
        valueFraction_*refValue_
      +
        (1.0 - valueFraction_)
       *transform(I - sqr(nHat), this->patchInternalField())
    );

    parent_bctype::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::partialSlipFvPatchField<Type>::snGradTransformDiag() const
{
    const vectorField nHat(this->patch().nf());
    vectorField diag(nHat.size());

    diag.replace(vector::X, mag(nHat.component(vector::X)));
    diag.replace(vector::Y, mag(nHat.component(vector::Y)));
    diag.replace(vector::Z, mag(nHat.component(vector::Z)));

    return
        valueFraction_*pTraits<Type>::one
      + (1.0 - valueFraction_)
       *transformFieldMask<Type>(pow<vector, pTraits<Type>::rank>(diag));
}


template<class Type>
void Foam::partialSlipFvPatchField<Type>::write(Ostream& os) const
{
    this->parent_bctype::write(os);
    refValue_.writeEntry("refValue", os);
    valueFraction_.writeEntry("valueFraction", os);

    if (writeValue_)
    {
        os.writeEntry("writeValue", "true");
        this->writeEntry("value", os);
    }
}


// ************************************************************************* //
