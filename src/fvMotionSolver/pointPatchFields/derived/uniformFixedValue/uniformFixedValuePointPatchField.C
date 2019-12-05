/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "uniformFixedValuePointPatchField.H"
#include "SubField.H"
#include "polyPatch.H"

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class Type>
const Foam::polyPatch&
Foam::uniformFixedValuePointPatchField<Type>::getPatch(const pointPatch& p)
{
    const polyMesh& mesh = p.boundaryMesh().mesh()();
    label patchi = mesh.boundaryMesh().findPatchID(p.name());

    if (patchi == -1)
    {
        FatalErrorInFunction
            << "Cannot use uniformFixedValue on patch " << p.name()
            << " since there is no underlying mesh patch" << exit(FatalError);
    }
    return mesh.boundaryMesh()[patchi];
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformFixedValuePointPatchField<Type>::
uniformFixedValuePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    fixedValuePointPatchField<Type>(p, iF),
    uniformValue_()
{}


template<class Type>
Foam::uniformFixedValuePointPatchField<Type>::
uniformFixedValuePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<Type>(p, iF, dict, false),
    uniformValue_
    (
        PatchFunction1<Type>::New
        (
            this->getPatch(p),
            "uniformValue",
            dict,
            false           // generate point values
        )
    )
{
    if (dict.found("value"))
    {
        fixedValuePointPatchField<Type>::operator==
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
Foam::uniformFixedValuePointPatchField<Type>::
uniformFixedValuePointPatchField
(
    const uniformFixedValuePointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<Type>(ptf, p, iF, mapper),
    uniformValue_(ptf.uniformValue_.clone(this->getPatch(p)))
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
Foam::uniformFixedValuePointPatchField<Type>::
uniformFixedValuePointPatchField
(
    const uniformFixedValuePointPatchField<Type>& ptf
)
:
    fixedValuePointPatchField<Type>(ptf),
    uniformValue_(ptf.uniformValue_.clone(this->getPatch(this->patch())))
{}


template<class Type>
Foam::uniformFixedValuePointPatchField<Type>::
uniformFixedValuePointPatchField
(
    const uniformFixedValuePointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    fixedValuePointPatchField<Type>(ptf, iF),
    uniformValue_(ptf.uniformValue_.clone(this->getPatch(this->patch())))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformFixedValuePointPatchField<Type>::autoMap
(
    const pointPatchFieldMapper& mapper
)
{
    fixedValuePointPatchField<Type>::autoMap(mapper);
    uniformValue_().autoMap(mapper);

    if (uniformValue_().constant())
    {
        // If mapper is not dependent on time we're ok to evaluate
        this->evaluate();
    }
}


template<class Type>
void Foam::uniformFixedValuePointPatchField<Type>::rmap
(
    const pointPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValuePointPatchField<Type>::rmap(ptf, addr);

    const uniformFixedValuePointPatchField& tiptf =
        refCast<const uniformFixedValuePointPatchField>(ptf);

    uniformValue_().rmap(tiptf.uniformValue_(), addr);
}


template<class Type>
void Foam::uniformFixedValuePointPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    const scalar t = this->db().time().timeOutputValue();
    fixedValuePointPatchField<Type>::operator==(uniformValue_->value(t));
    fixedValuePointPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::uniformFixedValuePointPatchField<Type>::
write(Ostream& os) const
{
    // Note: write value
    fixedValuePointPatchField<Type>::write(os);
    uniformValue_->writeData(os);
}


// ************************************************************************* //
