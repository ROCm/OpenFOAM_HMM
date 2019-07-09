/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2012-2016 OpenFOAM Foundation
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

#include "pressure.H"
#include "volFields.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(pressure, 0);
    addToRunTimeSelectionTable(functionObject, pressure, dictionary);
}
}

const Foam::Enum
<
    Foam::functionObjects::pressure::mode
>
Foam::functionObjects::pressure::modeNames
({
    { STATIC, "static" },
    { TOTAL, "total" },
    { ISENTROPIC, "isentropic" },
    { STATIC_COEFF, "staticCoeff" },
    { TOTAL_COEFF, "totalCoeff" },
});

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::word Foam::functionObjects::pressure::resultName() const
{
    word rName;

    if (mode_ & STATIC)
    {
        rName = "static(" + fieldName_ + ")";
    }
    else if (mode_ & TOTAL)
    {
        rName = "total(" + fieldName_ + ")";
    }
    else if (mode_ & ISENTROPIC)
    {
        rName = "isentropic(" + fieldName_ + ")";
    }
    else
    {
        FatalErrorInFunction
            << "Unhandled calculation mode " << modeNames[mode_]
            << abort(FatalError);
    }

    if (mode_ & COEFF)
    {
        rName += "_coeff";
    }

    return rName;
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::pressure::rhoScale
(
    const volScalarField& p
) const
{
    if (p.dimensions() == dimPressure)
    {
        return tmp<volScalarField>::New
        (
            IOobject
            (
                "rhoScale",
                p.mesh().time().timeName(),
                p.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            p,
            fvPatchField<scalar>::calculatedType()
        );
    }

    if (!rhoInfInitialised_)
    {
        FatalErrorInFunction
            << type() << " " << name() << ": "
            << "pressure identified as incompressible, but reference "
            << "density is not set.  Please set 'rho' to 'rhoInf', and "
            << "set an appropriate value for 'rhoInf'"
            << exit(FatalError);
    }

    return dimensionedScalar("rhoInf", dimDensity, rhoInf_)*p;
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::pressure::rhoScale
(
    const volScalarField& p,
    const tmp<volScalarField>& tsf
) const
{
    if (p.dimensions() == dimPressure)
    {
        return lookupObject<volScalarField>(rhoName_)*tsf;
    }

    return dimensionedScalar("rhoInf", dimDensity, rhoInf_)*tsf;
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::pressure::calcPressure
(
    const volScalarField& p,
    const tmp<volScalarField>& tp
) const
{
    if (mode_ & TOTAL)
    {
        return
            tp
          + dimensionedScalar("pRef", dimPressure, pRef_)
          + rhoScale(p, 0.5*magSqr(lookupObject<volVectorField>(UName_)));
    }

    if (mode_ & ISENTROPIC)
    {
        const basicThermo* thermoPtr =
            p.mesh().lookupObjectPtr<basicThermo>(basicThermo::dictName);

        if (!thermoPtr)
        {
            FatalErrorInFunction
                << "Isentropic pressure calculation requires a "
                << "thermodynamics package"
                << exit(FatalError);
        }

        const volScalarField gamma(thermoPtr->gamma());
        const volScalarField Mb
        (
            mag(lookupObject<volVectorField>(UName_))
           /sqrt(gamma*tp.ref()/thermoPtr->rho())
        );

        return tp()*(pow(1 + (gamma - 1)/2*sqr(Mb), gamma/(gamma - 1)));
    }

    return tp + dimensionedScalar("pRef", dimPressure, pRef_);
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::pressure::coeff
(
    const tmp<volScalarField>& tp
) const
{
    if (mode_ & COEFF)
    {
        tmp<volScalarField> tpCoeff(tp.ptr());
        volScalarField& pCoeff = tpCoeff.ref();

        pCoeff -= dimensionedScalar("pInf", dimPressure, pInf_);

        const dimensionedScalar pSmall("pSmall", dimPressure, SMALL);
        const dimensionedVector U("U", dimVelocity, UInf_);
        const dimensionedScalar rho("rho", dimDensity, rhoInf_);

        pCoeff /= 0.5*rho*magSqr(U) + pSmall;

        return tpCoeff;
    }

    return std::move(tp);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::pressure::calc()
{
    if (foundObject<volScalarField>(fieldName_))
    {
        const volScalarField& p = lookupObject<volScalarField>(fieldName_);

        auto tp = tmp<volScalarField>::New
        (
            IOobject
            (
                resultName_,
                p.mesh().time().timeName(),
                p.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            coeff(calcPressure(p, rhoScale(p)))
        );

        return store(resultName_, tp);
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::pressure::pressure
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "p"),
    mode_(STATIC),
    UName_("U"),
    rhoName_("rho"),
    pRef_(0),
    pInf_(0),
    UInf_(Zero),
    rhoInf_(1),
    rhoInfInitialised_(false)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::pressure::read(const dictionary& dict)
{
    Info<< type() << " " << name() << ":" << nl;

    fieldExpression::read(dict);

    UName_   = dict.lookupOrDefault<word>("U", "U");
    rhoName_ = dict.lookupOrDefault<word>("rho", "rho");

    if (rhoName_ == "rhoInf")
    {
        dict.readEntry("rhoInf", rhoInf_);
        rhoInfInitialised_ = true;
    }

    if (!modeNames.readIfPresent("mode", dict, mode_))
    {
        // Backwards compatibility
        // - check for the presence of 'calcTotal' and 'calcCoeff'

        bool calcTotal =
            dict.getOrDefaultCompat<bool>("mode", {{"calcTotal", 1812}}, false);
        bool calcCoeff =
            dict.getOrDefaultCompat<bool>("mode", {{"calcCoeff", 1812}}, false);

        if (calcTotal)
        {
            mode_ = TOTAL;
        }
        else
        {
            mode_ = STATIC;
        }

        if (calcCoeff)
        {
            mode_ = static_cast<mode>(COEFF | mode_);
        }
    }

    Info<< "    operating mode: " << modeNames[mode_] << nl;

    pRef_ = dict.lookupOrDefault<scalar>("pRef", 0);

    if (mode_ & COEFF)
    {
        dict.readEntry("pInf", pInf_);
        dict.readEntry("UInf", UInf_);
        dict.readEntry("rhoInf", rhoInf_);

        scalar zeroCheck = 0.5*rhoInf_*magSqr(UInf_) + pInf_;

        if (mag(zeroCheck) < ROOTVSMALL)
        {
            WarningInFunction
                << type() << " " << name() << ": "
                << "Coefficient calculation requested, but reference "
                << "pressure level is zero.  Please check the supplied "
                << "values of pInf, UInf and rhoInf" << endl;
        }

        rhoInfInitialised_ = true;
    }

    resultName_ = dict.lookupOrDefault<word>("result", resultName());

    Info<< endl;

    return true;
}


// ************************************************************************* //
