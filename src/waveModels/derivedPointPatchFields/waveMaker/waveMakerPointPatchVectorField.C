/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2018-2019 IH-Cantabria
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

#include "waveMakerPointPatchVectorField.H"
#include "mathematicalConstants.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "gravityMeshObject.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

const Foam::Enum<Foam::waveMakerPointPatchVectorField::motionTypes>
Foam::waveMakerPointPatchVectorField::motionTypeNames
({
    { motionTypes::piston, "piston" },
    { motionTypes::flap, "flap" },
    { motionTypes::solitary, "solitary" }
});


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const Foam::vector& Foam::waveMakerPointPatchVectorField::g()
{
    const meshObjects::gravity& gf = meshObjects::gravity::New(db().time());

    if (mag(gf.value()) < SMALL)
    {
        FatalErrorInFunction
            << "Gravity vector is not set.  Please update "
            << gf.uniformDimensionedVectorField::path()
            << exit(FatalError);
    }

    return gf.value();
}


Foam::scalar Foam::waveMakerPointPatchVectorField::waveLength
(
    const scalar h,
    const scalar T
)
{
    const scalar L0 = mag(g())*T*T/(constant::mathematical::twoPi);
    scalar L = L0;

    for (label i=1; i<=100; ++i)
    {
        L = L0*tanh(constant::mathematical::twoPi*h/L);
    }

    return L;
}


Foam::scalar Foam::waveMakerPointPatchVectorField::timeCoeff
(
    const scalar t
) const
{
    return max(0, min(t/rampTime_, 1));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveMakerPointPatchVectorField::waveMakerPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    motionType_(motionTypes::piston),
    n_(Zero),
    gHat_(Zero),
    initialDepth_(0),
    wavePeriod_(0),
    waveHeight_(0),
    wavePhase_(0),
    waveLength_(0),
    startTime_(0),
    rampTime_(1),
    secondOrder_(false)
{}


Foam::waveMakerPointPatchVectorField::waveMakerPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict, false),
    motionType_(motionTypeNames.get("motionType", dict)),
    n_(dict.get<vector>("n")),
    gHat_(Zero),
    initialDepth_(dict.get<scalar>("initialDepth")),
    wavePeriod_(dict.get<scalar>("wavePeriod")),
    waveHeight_(dict.get<scalar>("waveHeight")),
    wavePhase_(dict.get<scalar>("wavePhase")),
    waveLength_(this->waveLength(initialDepth_, wavePeriod_)),
    startTime_
    (
        dict.lookupOrDefault<scalar>
        (
            "startTime",
            db().time().startTime().value()
        )
    ),
    rampTime_(dict.get<scalar>("rampTime")),
    secondOrder_(dict.lookupOrDefault<bool>("secondOrder", false))
{
    // Create the co-ordinate system
    if (mag(n_) < SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Patch normal direction vector is not set.  'n' = " << n_
            << exit(FatalIOError);
    }
    n_ /= mag(n_);

    gHat_ = (g() - n_*(n_&g()));
    gHat_ /= mag(gHat_);

    if (!dict.found("value"))
    {
        updateCoeffs();
    }
}


Foam::waveMakerPointPatchVectorField::waveMakerPointPatchVectorField
(
    const waveMakerPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    motionType_(ptf.motionType_),
    n_(ptf.n_),
    gHat_(ptf.gHat_),
    initialDepth_(ptf.initialDepth_),
    wavePeriod_(ptf.wavePeriod_),
    waveHeight_(ptf.waveHeight_),
    wavePhase_(ptf.wavePhase_),
    waveLength_(ptf.waveLength_),
    startTime_(ptf.startTime_),
    rampTime_(ptf.rampTime_),
    secondOrder_(ptf.secondOrder_)
{}


Foam::waveMakerPointPatchVectorField::waveMakerPointPatchVectorField
(
    const waveMakerPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    motionType_(ptf.motionType_),
    n_(ptf.n_),
    gHat_(ptf.gHat_),
    initialDepth_(ptf.initialDepth_),
    wavePeriod_(ptf.wavePeriod_),
    waveHeight_(ptf.waveHeight_),
    wavePhase_(ptf.wavePhase_),
    waveLength_(ptf.waveLength_),
    startTime_(ptf.startTime_),
    rampTime_(ptf.rampTime_),
    secondOrder_(ptf.secondOrder_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::waveMakerPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalar t = db().time().value() - startTime_;

    const scalar waveK = constant::mathematical::twoPi/waveLength_;
    const scalar sigma = constant::mathematical::twoPi/wavePeriod_;

    const scalar kh = waveK*initialDepth_;

    switch (motionType_)
    {
        case motionTypes::flap:
        {
            const scalar m1 =
                4*sinh(kh)/(sinh(2*kh) + 2*kh)*(sinh(kh) + (1 - cosh(kh))/kh);

            scalar motionX = 0.5*waveHeight_/m1*sin(sigma*t);

            if (secondOrder_)
            {
                motionX +=
                    sqr(waveHeight_)/(16*initialDepth_)
                   *(3*cosh(kh)/pow3(sinh(kh)) - 2/m1)
                   *sin(2*sigma*t);
            }

            const pointField& points = patch().localPoints();
            const scalarField dz(-(points & gHat_) - initialDepth_);

            Field<vector>::operator=
            (
                n_*timeCoeff(t)*motionX*(1 + dz/initialDepth_)
            );

            break;
        }
        case motionTypes::piston:
        {
            const scalar m1 = 2*(cosh(2*kh) - 1)/(sinh(2*kh) + 2*kh);

            scalar motionX = 0.5*waveHeight_/m1*sin(sigma*t);

            if (secondOrder_)
            {
                motionX +=
                    sqr(waveHeight_)
                   /(32*initialDepth_)*(3*cosh(kh)
                   /pow3(sinh(kh)) - 2/m1);
            }

            Field<vector>::operator=(n_*timeCoeff(t)*motionX);

            break;
        }
        case motionTypes::solitary:
        {
            const scalar kappa = sqrt(0.75*waveHeight_/(pow3(initialDepth_)));
            const scalar waveCelerity =
                sqrt(mag(g())*(initialDepth_ + waveHeight_));
            const scalar stroke = sqrt(16.0*waveHeight_*initialDepth_/3.0);
            const scalar hr = waveHeight_/initialDepth_;
            wavePeriod_ = (2.0/(kappa*waveCelerity))*(3.8 + hr);
            const scalar tSolitary = -0.5*wavePeriod_ + t;

            // Newton-Rapshon
            scalar theta1 = 0;
            scalar theta2 = 0;
            scalar er = 10000;
            const scalar error = 0.001;
            while (er > error)
            {
                theta2 =
                    theta1
                  - (theta1 - kappa*waveCelerity*tSolitary + hr*tanh(theta1))
                   /(1.0 + hr*(1.0/cosh(theta1))*(1.0/cosh(theta1)));

                    er = mag(theta1 - theta2);
                    theta1 = theta2;
            }

            scalar motionX =
                waveHeight_/(kappa*initialDepth_)*tanh(theta1) + 0.5*stroke;

            Field<vector>::operator=(n_*motionX);

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled enumeration " << motionTypeNames[motionType_]
                << abort(FatalError);
        }
    }


    fixedValuePointPatchField<vector>::updateCoeffs();
}


void Foam::waveMakerPointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    os.writeEntry("motionType", motionTypeNames[motionType_]);
    os.writeEntry("n", n_);
    os.writeEntry("initialDepth", initialDepth_);
    os.writeEntry("wavePeriod", wavePeriod_);
    os.writeEntry("waveHeight", waveHeight_);
    os.writeEntry("wavePhase", wavePhase_);
    os.writeEntry("startTime", startTime_);
    os.writeEntry("rampTime", rampTime_);
    os.writeEntry("secondOrder", secondOrder_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePointPatchTypeField
    (
        pointPatchVectorField,
        waveMakerPointPatchVectorField
    );
}

// ************************************************************************* //
