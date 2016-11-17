/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
     \\/     M anipulation  | Copyright (C) 2015 IH-Cantabria
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

#include "waveInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "waveModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveInletVelocityFvPatchVectorField::waveInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    waveDictName_(waveModel::dictName)
{}


Foam::waveInletVelocityFvPatchVectorField::waveInletVelocityFvPatchVectorField
(
    const waveInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    waveDictName_(ptf.waveDictName_)
{}


Foam::waveInletVelocityFvPatchVectorField::waveInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    waveDictName_(dict.lookupOrDefault<word>("waveDict", waveModel::dictName))
{}


Foam::waveInletVelocityFvPatchVectorField::waveInletVelocityFvPatchVectorField
(
    const waveInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    waveDictName_(ptf.waveDictName_)
{}


Foam::waveInletVelocityFvPatchVectorField::waveInletVelocityFvPatchVectorField
(
    const waveInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    waveDictName_(ptf.waveDictName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::waveInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    tmp<waveModel> tmodel
    (
        waveModel::lookupOrCreate
        (
            patch().patch(),
            internalField().mesh(),
            waveDictName_
        )
    );

    waveModel& model = const_cast<waveModel&>(tmodel());

    model.correct(db().time().value());

    operator == (model.U());

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::waveInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);

    os.writeKeyword("waveDictName") << waveDictName_
        << token::END_STATEMENT << nl;

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       waveInletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
