/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016, 2019 OpenFOAM Foundation
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

#include "omegaWallFunctionFvPatchScalarField.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::omegaWallFunctionFvPatchScalarField::blendingType
>
Foam::omegaWallFunctionFvPatchScalarField::blendingTypeNames
({
    { blendingType::STEPWISE , "stepwise" },
    { blendingType::MAX , "max" },
    { blendingType::BINOMIAL2 , "binomial2" },
    { blendingType::BINOMIAL , "binomial" },
    { blendingType::EXPONENTIAL, "exponential" },
    { blendingType::TANH, "tanh" }
});

Foam::scalar Foam::omegaWallFunctionFvPatchScalarField::tolerance_ = 1e-5;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::omegaWallFunctionFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const volScalarField& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    label master = -1;
    forAll(bf, patchi)
    {
        if (isA<omegaWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            omegaWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            if (master == -1)
            {
                master = patchi;
            }

            opf.master() = master;
        }
    }
}


void Foam::omegaWallFunctionFvPatchScalarField::createAveragingWeights()
{
    const volScalarField& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    const fvMesh& mesh = omega.mesh();

    if (initialised_ && !mesh.changing())
    {
        return;
    }

    volScalarField weights
    (
        IOobject
        (
            "weights",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false // do not register
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    );

    DynamicList<label> omegaPatches(bf.size());
    forAll(bf, patchi)
    {
        if (isA<omegaWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            omegaPatches.append(patchi);

            const labelUList& faceCells = bf[patchi].patch().faceCells();
            for (const auto& celli : faceCells)
            {
                ++weights[celli];
            }
        }
    }

    cornerWeights_.setSize(bf.size());
    for (const auto& patchi : omegaPatches)
    {
        const fvPatchScalarField& wf = weights.boundaryField()[patchi];
        cornerWeights_[patchi] = 1.0/wf.patchInternalField();
    }

    G_.setSize(internalField().size(), 0.0);
    omega_.setSize(internalField().size(), 0.0);

    initialised_ = true;
}


Foam::omegaWallFunctionFvPatchScalarField&
Foam::omegaWallFunctionFvPatchScalarField::omegaPatch
(
    const label patchi
)
{
    const volScalarField& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    const omegaWallFunctionFvPatchScalarField& opf =
        refCast<const omegaWallFunctionFvPatchScalarField>(bf[patchi]);

    return const_cast<omegaWallFunctionFvPatchScalarField&>(opf);
}


void Foam::omegaWallFunctionFvPatchScalarField::calculateTurbulenceFields
(
    const turbulenceModel& turbModel,
    scalarField& G0,
    scalarField& omega0
)
{
    // accumulate all of the G and omega contributions
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            omegaWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            const List<scalar>& w = cornerWeights_[patchi];

            opf.calculate(turbModel, w, opf.patch(), G0, omega0);
        }
    }

    // apply zero-gradient condition for omega
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            omegaWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            opf == scalarField(omega0, opf.patch().faceCells());
        }
    }
}


void Foam::omegaWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel& turbModel,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& G0,
    scalarField& omega0
)
{
    const label patchi = patch.index();

    const nutWallFunctionFvPatchScalarField& nutw =
        nutWallFunctionFvPatchScalarField::nutw(turbModel, patchi);

    const scalarField& y = turbModel.y()[patchi];

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];

    const scalarField magGradUw(mag(Uw.snGrad()));

    const scalar Cmu25 = pow025(nutw.Cmu());

    // Set omega and G
    forAll(nutw, facei)
    {
        const label celli = patch.faceCells()[facei];
        const scalar yPlus = Cmu25*y[facei]*sqrt(k[celli])/nuw[facei];
        const scalar w = cornerWeights[facei];

        // Contribution from the viscous sublayer
        const scalar omegaVis = 6.0*nuw[facei]/(beta1_*sqr(y[facei]));

        // Contribution from the inertial sublayer
        const scalar omegaLog = sqrt(k[celli])/(Cmu25*nutw.kappa()*y[facei]);

        switch (blending_)
        {
            case blendingType::STEPWISE:
            {
                if (yPlus > nutw.yPlusLam())
                {
                    omega0[celli] += w*omegaLog;
                }
                else
                {
                    omega0[celli] += w*omegaVis;
                }
                break;
            }

            case blendingType::MAX:
            {
                // (PH:Eq. 27)
                omega0[celli] += max(omegaVis, omegaLog);
                break;
            }

            case blendingType::BINOMIAL2:
            {
                // (ME:Eq. 15)
                omega0[celli] += w*sqrt(sqr(omegaVis) + sqr(omegaLog));
                break;
            }

            case blendingType::BINOMIAL:
            {
                omega0[celli] +=
                    w*pow
                    (
                        pow(omegaVis, n_) + pow(omegaLog, n_),
                        1.0/n_
                    );
                break;
            }

            case blendingType::EXPONENTIAL:
            {
                // (PH:Eq. 31)
                const scalar Gamma = 0.01*pow4(yPlus)/(1.0 + 5.0*yPlus);
                const scalar invGamma = 1.0/(Gamma + ROOTVSMALL);

                omega0[celli] +=
                    w*(omegaVis*exp(-Gamma) + omegaLog*exp(-invGamma));
                break;
            }

            case blendingType::TANH:
            {
                // (KAS:Eqs. 33-34)
                const scalar phiTanh = tanh(pow4(yPlus/10.0));
                const scalar omegab1 = omegaVis + omegaLog;
                const scalar omegab2 =
                    pow(pow(omegaVis, 1.2) + pow(omegaLog, 1.2), 1.0/1.2);

                omega0[celli] += phiTanh*omegab1 + (1.0 - phiTanh)*omegab2;
                break;
            }
        }

        if (!(blending_ == blendingType::STEPWISE) || yPlus > nutw.yPlusLam())
        {
            G0[celli] +=
                w
               *(nutw[facei] + nuw[facei])
               *magGradUw[facei]
               *Cmu25*sqrt(k[celli])
               /(nutw.kappa()*y[facei]);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    blending_(blendingType::BINOMIAL2),
    n_(2.0),
    initialised_(false),
    master_(-1),
    beta1_(0.075),
    G_(),
    omega_(),
    cornerWeights_()
{}


Foam::omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const omegaWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    blending_(ptf.blending_),
    n_(ptf.n_),
    initialised_(false),
    master_(-1),
    beta1_(ptf.beta1_),
    G_(),
    omega_(),
    cornerWeights_()
{}


Foam::omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    blending_
    (
        blendingTypeNames.getOrDefault
        (
            "blending",
            dict,
            blendingType::BINOMIAL2
        )
    ),
    n_
    (
        dict.getCheckOrDefault<scalar>
        (
            "n",
            2.0,
            scalarMinMax::ge(0)
        )
    ),
    initialised_(false),
    master_(-1),
    beta1_(dict.getOrDefault<scalar>("beta1", 0.075)),
    G_(),
    omega_(),
    cornerWeights_()
{
    // The deprecated 'blended' keyword is superseded by the enum 'blending'
    if (dict.found("blended"))
    {
        IOWarningInFunction(dict)
            << "Using deprecated 'blended' keyword"
            << nl << "    Please use either of the below for the same behaviour:"
            << nl << "    'blending  binomial2;' for 'blended  on;'"
            << nl << "    'blending  stepwise;'  for 'blended  off;'"
            << nl << "    OVERWRITING: 'blended' keyword -> 'blending' enum"
            << endl;

        bool blended = dict.get<bool>("blended");

        if (blended)
        {
            blending_ = blendingType::BINOMIAL2;
        }
        else
        {
            blending_ = blendingType::STEPWISE;
        }
    }

    // apply zero-gradient condition on start-up
    this->operator==(patchInternalField());
}


Foam::omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const omegaWallFunctionFvPatchScalarField& owfpsf
)
:
    fixedValueFvPatchField<scalar>(owfpsf),
    blending_(owfpsf.blending_),
    n_(owfpsf.n_),
    initialised_(false),
    master_(-1),
    beta1_(owfpsf.beta1_),
    G_(),
    omega_(),
    cornerWeights_()
{}


Foam::omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const omegaWallFunctionFvPatchScalarField& owfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(owfpsf, iF),
    blending_(owfpsf.blending_),
    n_(owfpsf.n_),
    initialised_(false),
    master_(-1),
    beta1_(owfpsf.beta1_),
    G_(),
    omega_(),
    cornerWeights_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalarField& Foam::omegaWallFunctionFvPatchScalarField::G
(
    bool init
)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            G_ = 0.0;
        }

        return G_;
    }

    return omegaPatch(master_).G();
}


Foam::scalarField& Foam::omegaWallFunctionFvPatchScalarField::omega
(
    bool init
)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            omega_ = 0.0;
        }

        return omega_;
    }

    return omegaPatch(master_).omega(init);
}


void Foam::omegaWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbModel, G(true), omega(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& omega0 = this->omega();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G = db().lookupObjectRef<FieldType>(turbModel.GName());

    FieldType& omega = const_cast<FieldType&>(internalField());

    forAll(*this, facei)
    {
        const label celli = patch().faceCells()[facei];

        G[celli] = G0[celli];
        omega[celli] = omega0[celli];
    }

    fvPatchField<scalar>::updateCoeffs();
}


void Foam::omegaWallFunctionFvPatchScalarField::updateWeightedCoeffs
(
    const scalarField& weights
)
{
    if (updated())
    {
        return;
    }

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbModel, G(true), omega(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& omega0 = this->omega();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G = db().lookupObjectRef<FieldType>(turbModel.GName());

    FieldType& omega = const_cast<FieldType&>(internalField());

    scalarField& omegaf = *this;

    // only set the values if the weights are > tolerance
    forAll(weights, facei)
    {
        const scalar w = weights[facei];

        if (w > tolerance_)
        {
            const label celli = patch().faceCells()[facei];

            G[celli] = (1.0 - w)*G[celli] + w*G0[celli];
            omega[celli] = (1.0 - w)*omega[celli] + w*omega0[celli];
            omegaf[facei] = omega[celli];
        }
    }

    fvPatchField<scalar>::updateCoeffs();
}


void Foam::omegaWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    matrix.setValues(patch().faceCells(), patchInternalField());

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void Foam::omegaWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix,
    const Field<scalar>& weights
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    DynamicList<label> constraintCells(weights.size());
    DynamicList<scalar> constraintValues(weights.size());
    const labelUList& faceCells = patch().faceCells();

    const DimensionedField<scalar, volMesh>& fld = internalField();

    forAll(weights, facei)
    {
        // only set the values if the weights are > tolerance
        if (tolerance_ < weights[facei])
        {
            const label celli = faceCells[facei];

            constraintCells.append(celli);
            constraintValues.append(fld[celli]);
        }
    }

    if (debug)
    {
        Pout<< "Patch: " << patch().name()
            << ": number of constrained cells = " << constraintCells.size()
            << " out of " << patch().size()
            << endl;
    }

    matrix.setValues(constraintCells, constraintValues);

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void Foam::omegaWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    os.writeEntry("blending", blendingTypeNames[blending_]);
    os.writeEntry("n", n_);
    os.writeEntry("beta1", beta1_);
    fixedValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        omegaWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
