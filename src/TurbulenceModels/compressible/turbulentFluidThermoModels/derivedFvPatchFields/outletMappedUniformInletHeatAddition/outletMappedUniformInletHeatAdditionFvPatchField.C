/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

#include "outletMappedUniformInletHeatAdditionFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::outletMappedUniformInletHeatAdditionFvPatchField::
outletMappedUniformInletHeatAdditionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    outletPatchName_(),
    phiName_("phi"),
    Q_(0),
    minTempLimit_(0),
    maxTempLimit_(5000)
{}


Foam::outletMappedUniformInletHeatAdditionFvPatchField::
outletMappedUniformInletHeatAdditionFvPatchField
(
    const outletMappedUniformInletHeatAdditionFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    outletPatchName_(ptf.outletPatchName_),
    phiName_(ptf.phiName_),
    Q_(ptf.Q_),
    minTempLimit_(ptf.minTempLimit_),
    maxTempLimit_(ptf.maxTempLimit_)
{}


Foam::outletMappedUniformInletHeatAdditionFvPatchField::
outletMappedUniformInletHeatAdditionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    outletPatchName_(dict.lookup("outletPatch")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    Q_(readScalar(dict.lookup("Q"))),
    minTempLimit_(dict.lookupOrDefault<scalar>("minTempLimit", 0)),
    maxTempLimit_(dict.lookupOrDefault<scalar>("maxTempLimit", 5000))
{}



Foam::outletMappedUniformInletHeatAdditionFvPatchField::
outletMappedUniformInletHeatAdditionFvPatchField
(
    const outletMappedUniformInletHeatAdditionFvPatchField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    outletPatchName_(ptf.outletPatchName_),
    phiName_(ptf.phiName_),
    Q_(ptf.Q_),
    minTempLimit_(ptf.minTempLimit_),
    maxTempLimit_(ptf.maxTempLimit_)
{}



Foam::outletMappedUniformInletHeatAdditionFvPatchField::
outletMappedUniformInletHeatAdditionFvPatchField
(
    const outletMappedUniformInletHeatAdditionFvPatchField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    outletPatchName_(ptf.outletPatchName_),
    phiName_(ptf.phiName_),
    Q_(ptf.Q_),
    minTempLimit_(ptf.minTempLimit_),
    maxTempLimit_(ptf.maxTempLimit_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::outletMappedUniformInletHeatAdditionFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const GeometricField<scalar, fvPatchField, volMesh>& f
    (
        dynamic_cast<const GeometricField<scalar, fvPatchField, volMesh>&>
        (
            this->internalField()
        )
    );

    const fvPatch& p = this->patch();

    label outletPatchID =
        p.patch().boundaryMesh().findPatchID(outletPatchName_);

    if (outletPatchID < 0)
    {
        FatalErrorInFunction
            << "Unable to find outlet patch " << outletPatchName_
            << abort(FatalError);
    }

    const fvPatch& outletPatch = p.boundaryMesh()[outletPatchID];

    const fvPatchField<scalar>& outletPatchField =
        f.boundaryField()[outletPatchID];

    const surfaceScalarField& phi =
        this->db().lookupObject<surfaceScalarField>
        (
            phiName_
        );

    const scalarField& outletPatchPhi = phi.boundaryField()[outletPatchID];
    scalar sumOutletPatchPhi = gSum(outletPatchPhi);

    if (sumOutletPatchPhi > SMALL)
    {
        const basicThermo& thermo =
             this->db().lookupObject<basicThermo>(basicThermo::dictName);

        scalar averageOutletField =
            gSum(outletPatchPhi*outletPatchField)/sumOutletPatchPhi;

        const scalarField Cpf(thermo.Cp()().boundaryField()[outletPatchID]);

        scalar totalPhiCp = gSum(outletPatchPhi)*gAverage(Cpf);

        operator==
        (
            min
            (
                max
                (
                    averageOutletField + Q_/totalPhiCp,
                    minTempLimit_
                ),
                maxTempLimit_
            )
        );
    }
    else
    {
        scalar averageOutletField =
            gSum(outletPatch.magSf()*outletPatchField)
           /gSum(outletPatch.magSf());

        operator==(averageOutletField);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::outletMappedUniformInletHeatAdditionFvPatchField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("outletPatch")
        << outletPatchName_ << token::END_STATEMENT << nl;

    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);

    os.writeKeyword("Q") << Q_ << token::END_STATEMENT << nl;
    os.writeKeyword("minTempLimit")
        << minTempLimit_ << token::END_STATEMENT << nl;
    os.writeKeyword("maxTempLimit")
        << maxTempLimit_ << token::END_STATEMENT << nl;

    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        outletMappedUniformInletHeatAdditionFvPatchField
    );
}


// ************************************************************************* //
