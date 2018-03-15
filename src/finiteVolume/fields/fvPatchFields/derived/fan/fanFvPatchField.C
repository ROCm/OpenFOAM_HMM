/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fanFvPatchField<Type>::calcFanJump()
{
    if (this->cyclicPatch().owner())
    {
        this->jump_ = this->jumpTable_->value(this->db().time().value());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fanFvPatchField<Type>::fanFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    uniformJumpFvPatchField<Type>(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    uniformJump_(false),
    nonDimensional_(false),
    rpm_(0.0),
    dm_(0.0)
{}


template<class Type>
Foam::fanFvPatchField<Type>::fanFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    uniformJumpFvPatchField<Type>(p, iF, dict),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    uniformJump_(dict.lookupOrDefault<bool>("uniformJump", false)),
    nonDimensional_(dict.lookupOrDefault<Switch>("nonDimensional", false)),
    rpm_(dict.lookupOrDefault<scalar>("rpm", 0.0)),
    dm_(dict.lookupOrDefault<scalar>("dm", 0.0))
{
    if (nonDimensional_)
    {
        dict.lookup("rpm") >> rpm_;
        dict.lookup("dm") >> dm_;
    }
}


template<class Type>
Foam::fanFvPatchField<Type>::fanFvPatchField
(
    const fanFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    uniformJumpFvPatchField<Type>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    uniformJump_(ptf.uniformJump_),
    nonDimensional_(ptf.nonDimensional_),
    rpm_(ptf.rpm_),
    dm_(ptf.dm_)
{}


template<class Type>
Foam::fanFvPatchField<Type>::fanFvPatchField
(
    const fanFvPatchField<Type>& ptf
)
:
    uniformJumpFvPatchField<Type>(ptf),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    uniformJump_(ptf.uniformJump_),
    nonDimensional_(ptf.nonDimensional_),
    rpm_(ptf.rpm_),
    dm_(ptf.dm_)
{}


template<class Type>
Foam::fanFvPatchField<Type>::fanFvPatchField
(
    const fanFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    uniformJumpFvPatchField<Type>(ptf, iF),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    uniformJump_(ptf.uniformJump_),
    nonDimensional_(ptf.nonDimensional_),
    rpm_(ptf.rpm_),
    dm_(ptf.dm_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fanFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    calcFanJump();

    // Call fixedJump variant - uniformJump will overwrite the jump value
    fixedJumpFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::fanFvPatchField<Type>::write(Ostream& os) const
{
    uniformJumpFvPatchField<Type>::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    os.writeEntryIfDifferent<bool>("uniformJump", false, uniformJump_);
    os.writeEntryIfDifferent<bool>("nonDimensional", false, nonDimensional_);
    os.writeEntryIfDifferent<scalar>("rpm", false, rpm_);
    os.writeEntryIfDifferent<scalar>("dm", false, dm_);
}


// ************************************************************************* //
