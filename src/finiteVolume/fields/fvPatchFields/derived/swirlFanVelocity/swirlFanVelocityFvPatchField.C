/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "swirlFanVelocityFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::swirlFanVelocityFvPatchField::calcFanJump()
{
    if (this->cyclicPatch().owner())
    {
        const scalar rpm = rpm_->value(this->db().time().timeOutputValue());

        const surfaceScalarField& phi =
            db().lookupObject<surfaceScalarField>(phiName_);

        const fvPatchField<scalar>& pOwner =
            patch().lookupPatchField<volScalarField, scalar>(pName_);

        const label nbrIndex = this->cyclicPatch().neighbPatchID();

        const fvPatch& nbrPatch = patch().boundaryMesh()[nbrIndex];

        const fvPatchField<scalar>& pSlave =
            nbrPatch.lookupPatchField<volScalarField, scalar>(pName_);

        scalarField deltaP(mag(pOwner - pSlave));

        if (phi.dimensions() == dimMass/dimTime)
        {
            deltaP /=
                patch().lookupPatchField<volScalarField, scalar>(rhoName_);
        }

        const vector axisHat =
            gSum(patch().nf()*patch().magSf())/gSum(patch().magSf());

        vectorField tanDir
        (
            axisHat ^ (patch().Cf() - origin_)
        );

        tanDir /= (mag(tanDir) + SMALL);

        scalarField magTangU(patch().size(), Zero);

        if (useRealRadius_)
        {
            const vectorField& pCf = patch().Cf();

            forAll(pCf, i)
            {
                const scalar rMag = mag(pCf[i] - origin_);

                if (rMag > rInner_ && rMag < rOuter_)
                {
                    magTangU[i] =
                        deltaP[i]/rMag/fanEff_/rpmToRads(rpm);
                }
            }
        }
        else
        {
            if (rEff_ <= 0)
            {
                FatalErrorInFunction
                    << "Effective radius rEff should be specified in the "<< nl
                    << "dictionary." << nl
                    << exit(FatalError);
            }
            magTangU =
                deltaP/rEff_/fanEff_/rpmToRads(rpm);
        }

        // Calculate the tangential velocity
        const vectorField tangentialVelocity(magTangU*tanDir);

        this->setJump(tangentialVelocity);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swirlFanVelocityFvPatchField::swirlFanVelocityFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedJumpFvPatchField<vector>(p, iF),
    phiName_("phi"),
    pName_("p"),
    rhoName_("rho"),
    origin_(),
    rpm_(),
    rEff_(0.0),
    fanEff_(1.0),
    rInner_(0.0),
    rOuter_(0.0),
    useRealRadius_(false)
{}


Foam::swirlFanVelocityFvPatchField::swirlFanVelocityFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedJumpFvPatchField<vector>(p, iF, dict),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    pName_(dict.getOrDefault<word>("p", "p")),
    rhoName_(dict.getOrDefault<word>("rho", "rho")),
    origin_
    (
        dict.getOrDefault
        (
            "origin",
            returnReduce(patch().size(), maxOp<label>())
          ? gSum(patch().Cf()*patch().magSf())/gSum(patch().magSf())
          : Zero
        )
    ),
    rpm_
    (
        this->cyclicPatch().owner()
      ? Function1<scalar>::New("rpm", dict)
      : nullptr
    ),
    rEff_(dict.getOrDefault<scalar>("rEff", 0)),
    fanEff_(dict.getOrDefault<scalar>("fanEff", 1)),
    rInner_(dict.getOrDefault<scalar>("rInner", 0)),
    rOuter_(dict.getOrDefault<scalar>("rOuter", 0)),
    useRealRadius_(dict.getOrDefault("useRealRadius", false))
{}


Foam::swirlFanVelocityFvPatchField::swirlFanVelocityFvPatchField
(
    const swirlFanVelocityFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedJumpFvPatchField<vector>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    pName_(ptf.pName_),
    rhoName_(ptf.rhoName_),
    origin_(ptf.origin_),
    rpm_(ptf.rpm_.clone()),
    rEff_(ptf.rEff_),
    fanEff_(ptf.fanEff_),
    rInner_(ptf.rInner_),
    rOuter_(ptf.rOuter_),
    useRealRadius_(ptf.useRealRadius_)
{}


Foam::swirlFanVelocityFvPatchField::swirlFanVelocityFvPatchField
(
    const swirlFanVelocityFvPatchField& ptf
)
:
    fixedJumpFvPatchField<vector>(ptf),
    phiName_(ptf.phiName_),
    pName_(ptf.pName_),
    rhoName_(ptf.rhoName_),
    origin_(ptf.origin_),
    rpm_(ptf.rpm_.clone()),
    rEff_(ptf.rEff_),
    rInner_(ptf.rInner_),
    rOuter_(ptf.rOuter_),
    useRealRadius_(ptf.useRealRadius_)
{}


Foam::swirlFanVelocityFvPatchField::swirlFanVelocityFvPatchField
(
    const swirlFanVelocityFvPatchField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedJumpFvPatchField<vector>(ptf, iF),
    phiName_(ptf.phiName_),
    pName_(ptf.pName_),
    rhoName_(ptf.rhoName_),
    origin_(ptf.origin_),
    rpm_(ptf.rpm_.clone()),
    rEff_(ptf.rEff_),
    rInner_(ptf.rInner_),
    rOuter_(ptf.rOuter_),
    useRealRadius_(ptf.useRealRadius_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::swirlFanVelocityFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    calcFanJump();
}


void Foam::swirlFanVelocityFvPatchField::write(Ostream& os) const
{
    fixedJumpFvPatchField<vector>::write(os);

    if (this->cyclicPatch().owner())
    {
        os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
        os.writeEntryIfDifferent<word>("p", "p", pName_);
        os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
        os.writeEntry("origin", origin_);

        if (rpm_)
        {
            rpm_->writeData(os);
        }

        os.writeEntryIfDifferent<scalar>("rEff", 0.0, rEff_);
        os.writeEntryIfDifferent<bool>("useRealRadius", false, useRealRadius_);
        os.writeEntryIfDifferent<scalar>("rInner", 0.0, rInner_);
        os.writeEntryIfDifferent<scalar>("rOuter", 0.0, rOuter_);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       swirlFanVelocityFvPatchField
   );
}

// ************************************************************************* //
