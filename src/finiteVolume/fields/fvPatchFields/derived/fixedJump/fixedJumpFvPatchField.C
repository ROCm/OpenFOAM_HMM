/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "fixedJumpFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fixedJumpFvPatchField<Type>::fixedJumpFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    jumpCyclicFvPatchField<Type>(p, iF),
    jump_(this->size(), Zero),
    jump0_(this->size(), Zero),
    relaxFactor_(-1),
    timeIndex_(-1)
{}


template<class Type>
Foam::fixedJumpFvPatchField<Type>::fixedJumpFvPatchField
(
    const fixedJumpFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    jumpCyclicFvPatchField<Type>(ptf, p, iF, mapper),
    jump_(ptf.jump_, mapper),
    jump0_(ptf.jump0_, mapper),
    relaxFactor_(ptf.relaxFactor_),
    timeIndex_(ptf.timeIndex_)
{}


template<class Type>
Foam::fixedJumpFvPatchField<Type>::fixedJumpFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    jumpCyclicFvPatchField<Type>(p, iF, dict),
    jump_(p.size(), Zero),
    jump0_(p.size(), Zero),
    relaxFactor_(dict.getOrDefault<scalar>("relax", -1)),
    timeIndex_(this->db().time().timeIndex())
{
    if (this->cyclicPatch().owner())
    {
        jump_ = Field<Type>("jump", dict, p.size());

        if (dict.found("jump0"))
        {
            jump0_ = Field<Type>("jump0", dict, p.size());
        }
    }

    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        this->evaluate(Pstream::commsTypes::blocking);
    }
}


template<class Type>
Foam::fixedJumpFvPatchField<Type>::fixedJumpFvPatchField
(
    const fixedJumpFvPatchField<Type>& ptf
)
:
    jumpCyclicFvPatchField<Type>(ptf),
    jump_(ptf.jump_),
    jump0_(ptf.jump0_),
    relaxFactor_(ptf.relaxFactor_),
    timeIndex_(ptf.timeIndex_)
{}


template<class Type>
Foam::fixedJumpFvPatchField<Type>::fixedJumpFvPatchField
(
    const fixedJumpFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    jumpCyclicFvPatchField<Type>(ptf, iF),
    jump_(ptf.jump_),
    jump0_(ptf.jump0_),
    relaxFactor_(ptf.relaxFactor_),
    timeIndex_(ptf.timeIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fixedJumpFvPatchField<Type>::jump() const
{
    if (this->cyclicPatch().owner())
    {
        return jump_;
    }
    else
    {
        return refCast<const fixedJumpFvPatchField<Type>>
        (
            this->neighbourPatchField()
        ).jump();
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fixedJumpFvPatchField<Type>::jump0() const
{
    if (this->cyclicPatch().owner())
    {
        return jump0_;
    }
    else
    {
        return refCast<const fixedJumpFvPatchField<Type>>
        (
            this->neighbourPatchField()
        ).jump0();
    }
}


template<class Type>
Foam::scalar Foam::fixedJumpFvPatchField<Type>::relaxFactor() const
{
    return relaxFactor_;
}


template<class Type>
void Foam::fixedJumpFvPatchField<Type>::relax()
{
    if (!this->cyclicPatch().owner() || relaxFactor_ < 0)
    {
        return;
    }

    jump_ = relaxFactor_*jump_ + (1 - relaxFactor_)*jump0_;

    if (timeIndex_ != this->db().time().timeIndex())
    {
        jump0_ = jump_;

        timeIndex_ = this->db().time().timeIndex();
    }
}


template<class Type>
void Foam::fixedJumpFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    jumpCyclicFvPatchField<Type>::autoMap(m);
    jump_.autoMap(m);
    jump0_.autoMap(m);
}


template<class Type>
void Foam::fixedJumpFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    jumpCyclicFvPatchField<Type>::rmap(ptf, addr);

    const auto& fjptf = refCast<const fixedJumpFvPatchField<Type>>(ptf);
    jump_.rmap(fjptf.jump_, addr);
    jump0_.rmap(fjptf.jump0_, addr);
}


template<class Type>
void Foam::fixedJumpFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    // Write patchType if not done already by fvPatchField
    if (!this->patchType().size())
    {
        os.writeEntry("patchType", this->interfaceFieldType());
    }

    if (this->cyclicPatch().owner())
    {
        jump_.writeEntry("jump", os);

        if (relaxFactor_ > 0)
        {
            os.writeEntry("relax", relaxFactor_);
            jump0_.writeEntry("jump0", os);
        }
    }

    this->writeEntry("value", os);
}


// ************************************************************************* //
