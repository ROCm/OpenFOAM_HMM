/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "fanFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fanFvPatchField<Type>::fanFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedJumpFvPatchField<Type>(p, iF),
    jumpTable_(0)
{}


template<class Type>
Foam::fanFvPatchField<Type>::fanFvPatchField
(
    const fanFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedJumpFvPatchField<Type>(ptf, p, iF, mapper),
    jumpTable_(ptf.jumpTable_().clone().ptr())
{}


template<class Type>
Foam::fanFvPatchField<Type>::fanFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedJumpFvPatchField<Type>(p, iF),
    jumpTable_(DataEntry<Type>::New("jumpTable", dict))
{}


template<class Type>
Foam::fanFvPatchField<Type>::fanFvPatchField
(
    const fanFvPatchField<Type>& ptf
)
:
    cyclicLduInterfaceField(),
    fixedJumpFvPatchField<Type>(ptf),
    jumpTable_(ptf.jumpTable_().clone().ptr())
{}


template<class Type>
Foam::fanFvPatchField<Type>::fanFvPatchField
(
    const fanFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedJumpFvPatchField<Type>(ptf, iF),
    jumpTable_(ptf.jumpTable_().clone().ptr())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
void Foam::fanFvPatchField<Type>::write(Ostream& os) const
{

    fixedJumpFvPatchField<Type>::write(os);
    if (this->cyclicPatch().owner())
    {
        jumpTable_->writeData(os);
    }
    this->writeEntry("value", os);
}


// ************************************************************************* //
