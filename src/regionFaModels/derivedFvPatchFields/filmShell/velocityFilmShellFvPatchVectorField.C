/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "velocityFilmShellFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

velocityFilmShellFvPatchVectorField::velocityFilmShellFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchField<vector>(p, iF),
    baffle_(),
    dict_(dictionary::null),
    curTimeIndex_(-1),
    zeroWallVelocity_(true)
{
    refValue() = 0;
    refGrad() = 0;
    valueFraction() = 1;
}


velocityFilmShellFvPatchVectorField::velocityFilmShellFvPatchVectorField
(
    const velocityFilmShellFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<vector>
    (
        ptf,
        p,
        iF,
        mapper
    ),
    baffle_(),
    dict_(ptf.dict_),
    curTimeIndex_(-1),
    zeroWallVelocity_(true)
{}


velocityFilmShellFvPatchVectorField::velocityFilmShellFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<vector>(p, iF),
    baffle_(nullptr),
    dict_(dict),
    curTimeIndex_(-1),
    zeroWallVelocity_(dict.getOrDefault<bool>("zeroWallVelocity", true))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

    typedef regionModels::areaSurfaceFilmModels::liquidFilmBase baffle;

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = vectorField("refValue", dict, p.size());
        refGrad() = vectorField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = vector::zero;
        valueFraction() = 1;
    }

    if (!baffle_)
    {
        baffle_.reset(baffle::New(p, dict).ptr());
    }
}


velocityFilmShellFvPatchVectorField::velocityFilmShellFvPatchVectorField
(
    const velocityFilmShellFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchField<vector>(ptf, iF),
    baffle_(),
    dict_(ptf.dict_),
    curTimeIndex_(-1),
    zeroWallVelocity_(true)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void velocityFilmShellFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Execute the change to the openFraction only once per time-step
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        baffle_->evolve();

        volVectorField::Boundary& vfb =
            db().lookupObjectRef<volVectorField>
            (
                this->internalField().name()
            ).boundaryFieldRef();

        baffle_->vsm().mapToVolume(baffle_->Us(), vfb);

        refGrad() = Zero;
        valueFraction() = 1;

        if (zeroWallVelocity_)
        {
            refValue() = Zero;
        }
        else
        {
            refValue() = vfb[patch().index()];
        }
        curTimeIndex_ = this->db().time().timeIndex();
    }

    mixedFvPatchField<vector>::updateCoeffs();
}


void velocityFilmShellFvPatchVectorField::write(Ostream& os) const
{
    mixedFvPatchField<vector>::write(os);

    // Remove value and type already written by mixedFvPatchField
    dict_.remove("value");
    dict_.remove("type");
    dict_.remove("refValue");
    dict_.remove("refGradient");
    dict_.remove("valueFraction");
    dict_.write(os, false);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    velocityFilmShellFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
