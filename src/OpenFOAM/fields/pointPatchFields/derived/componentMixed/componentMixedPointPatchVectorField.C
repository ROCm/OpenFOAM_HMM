/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "ComponentMixedPointPatchVectorField.H"
#include "PointPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class PointPatch
>
void ComponentMixedPointPatchVectorField<PatchField, PointPatch>
::checkFieldSize() const
{
    if
    (
        this->size() != this->patch().size()
     || refValue_.size() != this->patch().size()
     || valueFraction_.size() != this->patch().size()
    )
    {
        FatalErrorIn
        (
            "void ComponentMixedPointPatchVectorField::checkField() const"
        )   << "field does not correspond to patch. " << endl
            << "Field size: " << this->size()
            << " value size: " << refValue_.size()
            << " valueFraction size: " << valueFraction_.size()
            << " patch size: " << this->patch().size()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class PointPatch
>
ComponentMixedPointPatchVectorField<PatchField, PointPatch>::
ComponentMixedPointPatchVectorField
(
    const PointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    PatchField<vector>(p, iF),
    refValue_(p.size()),
    valueFraction_(p.size())
{}


template
<
    template<class> class PatchField,
    class PointPatch
>
ComponentMixedPointPatchVectorField<PatchField, PointPatch>::
ComponentMixedPointPatchVectorField
(
    const PointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    PatchField<vector>(p, iF, dict),
    refValue_("refValue", dict, p.size()),
    valueFraction_("valueFraction", dict, p.size())
{}


template
<
    template<class> class PatchField,
    class PointPatch
>
ComponentMixedPointPatchVectorField<PatchField, PointPatch>::
ComponentMixedPointPatchVectorField
(
    const ComponentMixedPointPatchVectorField<PatchField, PointPatch>& ptf,
    const PointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const PointPatchFieldMapper& mapper
)
:
    PatchField<vector>(ptf, iF),
    refValue_(ptf.refValue_, mapper),
    valueFraction_(ptf.valueFraction_, mapper)
{}


template
<
    template<class> class PatchField,
    class PointPatch
>
ComponentMixedPointPatchVectorField<PatchField, PointPatch>::
ComponentMixedPointPatchVectorField
(
    const ComponentMixedPointPatchVectorField<PatchField, PointPatch>& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    PatchField<vector>(ptf, iF),
    refValue_(ptf.refValue_),
    valueFraction_(ptf.valueFraction_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    template<class> class PatchField,
    class PointPatch
>
void ComponentMixedPointPatchVectorField<PatchField, PointPatch>::autoMap
(
    const PointPatchFieldMapper& m
)
{
    refValue_.autoMap(m);
    valueFraction_.autoMap(m);
}


// Grab the values using rmap
template
<
    template<class> class PatchField,
    class PointPatch
>
void ComponentMixedPointPatchVectorField<PatchField, PointPatch>::rmap
(
    const PointPatchField<PatchField, PointPatch, vector>& ptf,
    const labelList& addr
)
{
    const ComponentMixedPointPatchVectorField& mptf =
        refCast<const ComponentMixedPointPatchVectorField>(ptf);

    refValue_.rmap(mptf.refValue_, addr);
    valueFraction_.rmap(mptf.valueFraction_, addr);
}


// Evaluate patch field
template
<
    template<class> class PatchField,
    class PointPatch
>
void ComponentMixedPointPatchVectorField<PatchField, PointPatch>::evaluate
(
    const Pstream::commsTypes
)
{
    tmp<vectorField> internalValues = this->patchInternalField();

    // Get internal field to insert values into
    Field<vector>& iF = const_cast<vectorField&>(this->internalField());

    vectorField values =
        cmptMultiply(refValue_, valueFraction_)
      + cmptMultiply(internalValues, vector::one - valueFraction_);

    this->setInInternalField(iF, values);
}


// Write
template
<
    template<class> class PatchField,
    class PointPatch
>
void ComponentMixedPointPatchVectorField<PatchField, PointPatch>::
write(Ostream& os) const
{
    PatchField<vector>::write(os);
    refValue_.writeEntry("refValue", os);
    valueFraction_.writeEntry("valueFraction", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
