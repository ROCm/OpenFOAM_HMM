/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2017 OpenCFD Ltd.
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

#include "externalWallHeatFluxTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "physicoChemicalConstants.H"

using Foam::constant::physicoChemical::sigma;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::externalWallHeatFluxTemperatureFvPatchScalarField::operationMode
>
Foam::externalWallHeatFluxTemperatureFvPatchScalarField::operationModeNames
{
    { operationMode::fixedPower, "power" },
    { operationMode::fixedHeatFlux, "flux" },
    { operationMode::fixedHeatTransferCoeff, "coefficient" },
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalWallHeatFluxTemperatureFvPatchScalarField::
externalWallHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), "undefined", "undefined", "undefined-K"),
    mode_(fixedHeatFlux),
    Q_(0),
    Ta_(),
    relaxation_(1),
    emissivity_(0),
    qrRelaxation_(1),
    qrName_("undefined-qr"),
    thicknessLayers_(),
    kappaLayers_()
{
    refValue() = 0;
    refGrad() = 0;
    valueFraction() = 1;
}


Foam::externalWallHeatFluxTemperatureFvPatchScalarField::
externalWallHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    mode_(operationModeNames.lookup("mode", dict)),
    Q_(0),
    Ta_(Function1<scalar>::New("Ta", dict)),
    relaxation_(dict.lookupOrDefault<scalar>("relaxation", 1)),
    emissivity_(dict.lookupOrDefault<scalar>("emissivity", 0)),
    qrRelaxation_(dict.lookupOrDefault<scalar>("qrRelaxation", 1)),
    qrName_(dict.lookupOrDefault<word>("qr", "none")),
    thicknessLayers_(),
    kappaLayers_()
{
    switch (mode_)
    {
        case fixedPower:
        {
            dict.lookup("Q") >> Q_;

            break;
        }
        case fixedHeatFlux:
        {
            q_ = scalarField("q", dict, p.size());

            break;
        }
        case fixedHeatTransferCoeff:
        {
            h_ = scalarField("h", dict, p.size());

            if (dict.found("thicknessLayers"))
            {
                dict.lookup("thicknessLayers") >> thicknessLayers_;
                dict.lookup("kappaLayers") >> kappaLayers_;

                if (thicknessLayers_.size() != kappaLayers_.size())
                {
                    FatalIOErrorInFunction(dict)
                        << "\n number of layers for thicknessLayers and "
                        << "kappaLayers must be the same"
                        << "\n for patch " << p.name()
                        << " of field " << internalField().name()
                        << " in file " << internalField().objectPath()
                        << exit(FatalIOError);
                }
            }

            break;
        }
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (qrName_ != "none")
    {
        if (dict.found("qrPrevious"))
        {
            qrPrevious_ = scalarField("qrPrevious", dict, p.size());
        }
        else
        {
            qrPrevious_.setSize(p.size(), 0);
        }
    }

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0;
        valueFraction() = 1;
    }
}


Foam::externalWallHeatFluxTemperatureFvPatchScalarField::
externalWallHeatFluxTemperatureFvPatchScalarField
(
    const externalWallHeatFluxTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    temperatureCoupledBase(patch(), ptf),
    mode_(ptf.mode_),
    Q_(ptf.Q_),
    q_(ptf.q_, mapper),
    h_(ptf.h_, mapper),
    Ta_(ptf.Ta_, false),
    relaxation_(ptf.relaxation_),
    emissivity_(ptf.emissivity_),
    qrPrevious_(ptf.qrPrevious_, mapper),
    qrRelaxation_(ptf.qrRelaxation_),
    qrName_(ptf.qrName_),
    thicknessLayers_(ptf.thicknessLayers_),
    kappaLayers_(ptf.kappaLayers_)
{
    switch (mode_)
    {
        case fixedPower:
        {
            break;
        }
        case fixedHeatFlux:
        {
            q_.map(ptf.q_, mapper);

            break;
        }
        case fixedHeatTransferCoeff:
        {
            h_.map(ptf.h_, mapper);

            break;
        }
    }

    if (qrName_ != "none")
    {
        qrPrevious_.map(ptf.qrPrevious_, mapper);
    }
}


Foam::externalWallHeatFluxTemperatureFvPatchScalarField::
externalWallHeatFluxTemperatureFvPatchScalarField
(
    const externalWallHeatFluxTemperatureFvPatchScalarField& ewhftpsf
)
:
    mixedFvPatchScalarField(ewhftpsf),
    temperatureCoupledBase(ewhftpsf),
    mode_(ewhftpsf.mode_),
    Q_(ewhftpsf.Q_),
    q_(ewhftpsf.q_),
    h_(ewhftpsf.h_),
    Ta_(ewhftpsf.Ta_, false),
    relaxation_(ewhftpsf.relaxation_),
    emissivity_(ewhftpsf.emissivity_),
    qrPrevious_(ewhftpsf.qrPrevious_),
    qrRelaxation_(ewhftpsf.qrRelaxation_),
    qrName_(ewhftpsf.qrName_),
    thicknessLayers_(ewhftpsf.thicknessLayers_),
    kappaLayers_(ewhftpsf.kappaLayers_)
{}


Foam::externalWallHeatFluxTemperatureFvPatchScalarField::
externalWallHeatFluxTemperatureFvPatchScalarField
(
    const externalWallHeatFluxTemperatureFvPatchScalarField& ewhftpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ewhftpsf, iF),
    temperatureCoupledBase(patch(), ewhftpsf),
    mode_(ewhftpsf.mode_),
    Q_(ewhftpsf.Q_),
    q_(ewhftpsf.q_),
    h_(ewhftpsf.h_),
    Ta_(ewhftpsf.Ta_, false),
    relaxation_(ewhftpsf.relaxation_),
    emissivity_(ewhftpsf.emissivity_),
    qrPrevious_(ewhftpsf.qrPrevious_),
    qrRelaxation_(ewhftpsf.qrRelaxation_),
    qrName_(ewhftpsf.qrName_),
    thicknessLayers_(ewhftpsf.thicknessLayers_),
    kappaLayers_(ewhftpsf.kappaLayers_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::externalWallHeatFluxTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);

    switch (mode_)
    {
        case fixedPower:
        {
            break;
        }
        case fixedHeatFlux:
        {
            q_.autoMap(m);

            break;
        }
        case fixedHeatTransferCoeff:
        {
            h_.autoMap(m);

            break;
        }
    }

    if (qrName_ != "none")
    {
        qrPrevious_.autoMap(m);
    }
}


void Foam::externalWallHeatFluxTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const externalWallHeatFluxTemperatureFvPatchScalarField& ewhftpsf =
        refCast<const externalWallHeatFluxTemperatureFvPatchScalarField>(ptf);

    switch (mode_)
    {
        case fixedPower:
        {
            break;
        }
        case fixedHeatFlux:
        {
            q_.rmap(ewhftpsf.q_, addr);

            break;
        }
        case fixedHeatTransferCoeff:
        {
            h_.rmap(ewhftpsf.h_, addr);

            break;
        }
    }

    if (qrName_ != "none")
    {
        qrPrevious_.rmap(ewhftpsf.qrPrevious_, addr);
    }
}


void Foam::externalWallHeatFluxTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField& Tp(*this);

    scalarField qr(Tp.size(), 0);
    if (qrName_ != "none")
    {
        qr =
            qrRelaxation_
           *patch().lookupPatchField<volScalarField, scalar>(qrName_)
          + (1 - qrRelaxation_)*qrPrevious_;

        qrPrevious_ = qr;
    }

    switch (mode_)
    {
        case fixedPower:
        {
            refGrad() = (Q_/gSum(patch().magSf()) + qr)/kappa(Tp);
            refValue() = 0;
            valueFraction() = 0;

            break;
        }
        case fixedHeatFlux:
        {
            refGrad() = (q_ + qr)/kappa(Tp);
            refValue() = 0;
            valueFraction() = 0;

            break;
        }
        case fixedHeatTransferCoeff:
        {
            scalar totalSolidRes = 0;
            if (thicknessLayers_.size())
            {
                forAll(thicknessLayers_, iLayer)
                {
                    const scalar l = thicknessLayers_[iLayer];
                    if (kappaLayers_[iLayer] > 0)
                    {
                        totalSolidRes += l/kappaLayers_[iLayer];
                    }
                }
            }
            scalarField hp(1/(1/h_ + totalSolidRes));

            const scalar Ta = Ta_->value(this->db().time().timeOutputValue());
            scalarField hpTa(hp*Ta);

            if (emissivity_ > 0)
            {
                // Evaluate the radiative flux to the environment
                // from the surface temperature ...
                if (totalSolidRes > 0)
                {
                    // ... including the effect of the solid wall thermal
                    // resistance
                    scalarField TpLambda(h_/(h_ + 1/totalSolidRes));
                    scalarField Ts(TpLambda*Tp + (1 - TpLambda)*Ta);
                    scalarField lambdaTa4(pow4((1 - TpLambda)*Ta));

                    hp += emissivity_*sigma.value()*(pow4(Ts) - lambdaTa4)/Tp;
                    hpTa += sigma.value()*(emissivity_*lambdaTa4 + pow4(Ta));
                }
                else
                {
                    // ... if there is no solid wall thermal resistance use
                    // the current wall temperature
                    hp += emissivity_*sigma.value()*pow3(Tp);
                    hpTa += sigma.value()*pow4(Ta);
                }
            }

            const scalarField kappaDeltaCoeffs
            (
                this->kappa(Tp)*patch().deltaCoeffs()
            );

            refGrad() = 0;

            forAll(Tp, i)
            {
                if (qr[i] < 0)
                {
                    const scalar hpmqr = hp[i] - qr[i]/Tp[i];

                    refValue()[i] = hpTa[i]/hpmqr;
                    valueFraction()[i] = hpmqr/(hpmqr + kappaDeltaCoeffs[i]);
                }
                else
                {
                    refValue()[i] = (hpTa[i] + qr[i])/hp[i];
                    valueFraction()[i] = hp[i]/(hp[i] + kappaDeltaCoeffs[i]);
                }
            }

            break;
        }
    }

    valueFraction() = relaxation_*valueFraction() + (1 - relaxation_);
    refValue() = relaxation_*refValue() + (1 - relaxation_)*Tp;

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        const scalar Q = gSum(kappa(Tp)*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << internalField().name() << " :"
            << " heat transfer rate:" << Q
            << " wall temperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }
}


void Foam::externalWallHeatFluxTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);

    os.writeEntry("mode", operationModeNames[mode_]);
    temperatureCoupledBase::write(os);

    switch (mode_)
    {
        case fixedPower:
        {
            os.writeEntry("Q", Q_);
            break;
        }
        case fixedHeatFlux:
        {
            q_.writeEntry("q", os);

            break;
        }
        case fixedHeatTransferCoeff:
        {
            h_.writeEntry("h", os);
            Ta_->writeData(os);

            if (relaxation_ < 1)
            {
                os.writeEntry("relaxation", relaxation_);
            }

            if (emissivity_ > 0)
            {
                os.writeEntry("emissivity", emissivity_);
            }

            if (thicknessLayers_.size())
            {
                thicknessLayers_.writeEntry("thicknessLayers", os);
                kappaLayers_.writeEntry("kappaLayers", os);
            }

            break;
        }
    }

    os.writeEntry("qr", qrName_);

    if (qrName_ != "none")
    {
        os.writeEntry("qrRelaxation", qrRelaxation_);

        qrPrevious_.writeEntry("qrPrevious", os);
    }

    refValue().writeEntry("refValue", os);
    refGrad().writeEntry("refGradient", os);
    valueFraction().writeEntry("valueFraction", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        externalWallHeatFluxTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
