/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "turbulentTemperatureRadCoupledMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentTemperatureRadCoupledMixedFvPatchScalarField::
turbulentTemperatureRadCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase
    (
        patch(),
        "undefined",
        "undefined",
        "undefined-K",
        "undefined-alpha"
    ),
    TnbrName_("undefined-Tnbr"),
    qrNbrName_("undefined-qrNbr"),
    qrName_("undefined-qr"),
    thicknessLayers_(0),
    kappaLayers_(0),
    contactRes_(0),
    thermalInertia_(false)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


turbulentTemperatureRadCoupledMixedFvPatchScalarField::
turbulentTemperatureRadCoupledMixedFvPatchScalarField
(
    const turbulentTemperatureRadCoupledMixedFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    temperatureCoupledBase(patch(), psf),
    TnbrName_(psf.TnbrName_),
    qrNbrName_(psf.qrNbrName_),
    qrName_(psf.qrName_),
    thicknessLayers_(psf.thicknessLayers_),
    kappaLayers_(psf.kappaLayers_),
    contactRes_(psf.contactRes_),
    thermalInertia_(psf.thermalInertia_)
{}


turbulentTemperatureRadCoupledMixedFvPatchScalarField::
turbulentTemperatureRadCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    TnbrName_(dict.getOrDefault<word>("Tnbr", "T")),
    qrNbrName_(dict.getOrDefault<word>("qrNbr", "none")),
    qrName_(dict.getOrDefault<word>("qr", "none")),
    thicknessLayers_(0),
    kappaLayers_(0),
    contactRes_(0.0),
    thermalInertia_(dict.getOrDefault<Switch>("thermalInertia", false))
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

    if (dict.readIfPresent("thicknessLayers", thicknessLayers_))
    {
        dict.readEntry("kappaLayers", kappaLayers_);

        if (thicknessLayers_.size() > 0)
        {
            // Calculate effective thermal resistance by harmonic averaging
            forAll(thicknessLayers_, iLayer)
            {
                contactRes_ += thicknessLayers_[iLayer]/kappaLayers_[iLayer];
            }
            contactRes_ = 1.0/contactRes_;
        }
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

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
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


turbulentTemperatureRadCoupledMixedFvPatchScalarField::
turbulentTemperatureRadCoupledMixedFvPatchScalarField
(
    const turbulentTemperatureRadCoupledMixedFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    temperatureCoupledBase(patch(), psf),
    TnbrName_(psf.TnbrName_),
    qrNbrName_(psf.qrNbrName_),
    qrName_(psf.qrName_),
    thicknessLayers_(psf.thicknessLayers_),
    kappaLayers_(psf.kappaLayers_),
    contactRes_(psf.contactRes_),
    thermalInertia_(psf.thermalInertia_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentTemperatureRadCoupledMixedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const polyMesh& mesh = patch().boundaryMesh().mesh();

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the coupling information from the mappedPatchBase
    const label patchi = patch().index();
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const label samplePatchi = mpp.samplePolyPatch().index();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];


    scalarField Tc(patchInternalField());
    scalarField& Tp = *this;

    const turbulentTemperatureRadCoupledMixedFvPatchScalarField&
        nbrField = refCast
            <const turbulentTemperatureRadCoupledMixedFvPatchScalarField>
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>(TnbrName_)
            );

    // Swap to obtain full local values of neighbour K*delta
    scalarField KDeltaNbr;
    tmp<scalarField> TcNbr(new scalarField(nbrField.size(), Zero));
    if (contactRes_ == 0.0)
    {
        TcNbr.ref() = nbrField.patchInternalField();
        KDeltaNbr = nbrField.kappa(nbrField)*nbrPatch.deltaCoeffs();
    }
    else
    {
        TcNbr.ref() = nbrField;
        KDeltaNbr.setSize(nbrField.size(), contactRes_);
    }
    mpp.distribute(KDeltaNbr);
    mpp.distribute(TcNbr.ref());

    scalarField KDelta(kappa(Tp)*patch().deltaCoeffs());

    scalarField qr(Tp.size(), Zero);
    if (qrName_ != "none")
    {
        qr = patch().lookupPatchField<volScalarField, scalar>(qrName_);
    }

    scalarField qrNbr(Tp.size(), Zero);
    if (qrNbrName_ != "none")
    {
        qrNbr = nbrPatch.lookupPatchField<volScalarField, scalar>(qrNbrName_);
        mpp.distribute(qrNbr);
    }

    // inertia therm
    if (thermalInertia_)
    {
        const scalar dt = mesh.time().deltaTValue();
        scalarField mCpDtNbr;

        {
            const basicThermo* thermo =
                nbrMesh.findObject<basicThermo>(basicThermo::dictName);

            if (thermo)
            {
                const scalarField& ppn =
                    thermo->p().boundaryField()[samplePatchi];
                const scalarField& Tpn =
                    thermo->T().boundaryField()[samplePatchi];

                mCpDtNbr =
                (
                    thermo->Cp(ppn, Tpn, samplePatchi)
                  * thermo->rho(samplePatchi)
                  / nbrPatch.deltaCoeffs()/dt
                );

                mpp.distribute(mCpDtNbr);
            }
            else
            {
                mCpDtNbr.setSize(Tp.size(), Zero);
            }
        }

        scalarField mCpDt;

        // Local inertia therm
        {
            const basicThermo* thermo =
                mesh.findObject<basicThermo>(basicThermo::dictName);

            if (thermo)
            {
                const scalarField& pp = thermo->p().boundaryField()[patchi];

                mCpDt =
                (
                    thermo->Cp(pp, Tp, patchi)
                  * thermo->rho(patchi)
                  / patch().deltaCoeffs()/dt
                );
            }
            else
            {
                // Issue warning?
                mCpDt.setSize(Tp.size(), Zero);
            }
        }

        const volScalarField& T =
            this->db().lookupObject<volScalarField>
            (
                this->internalField().name()
            );

        const fvPatchField<scalar>& TpOld = T.oldTime().boundaryField()[patchi];

        scalarField alpha(KDeltaNbr + mCpDt + mCpDtNbr);

        valueFraction() = alpha/(alpha + KDelta);
        scalarField c(KDeltaNbr*TcNbr + (mCpDt + mCpDtNbr)*TpOld);
        refValue() = c/alpha;
        refGrad() = (qr + qrNbr)/kappa(Tp);
    }
    else
    {
        valueFraction() = KDeltaNbr/(KDeltaNbr + KDelta);
        refValue() = TcNbr;
        refGrad() = (qr + qrNbr)/kappa(Tp);
    }

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar Q = gSum(kappa(Tp)*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " <- "
            << nbrMesh.name() << ':'
            << nbrPatch.name() << ':'
            << this->internalField().name() << " :"
            << " heat transfer rate:" << Q
            << " walltemperature "
            << " min:" << gMin(Tp)
            << " max:" << gMax(Tp)
            << " avg:" << gAverage(Tp)
            << endl;
    }

    // Restore tag
    UPstream::msgType() = oldTag;
}


void turbulentTemperatureRadCoupledMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeEntry("Tnbr", TnbrName_);

    os.writeEntry("qrNbr", qrNbrName_);
    os.writeEntry("qr", qrName_);
    os.writeEntry("thermalInertia", thermalInertia_);

    thicknessLayers_.writeEntry("thicknessLayers", os);
    kappaLayers_.writeEntry("kappaLayers", os);

    temperatureCoupledBase::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentTemperatureRadCoupledMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
