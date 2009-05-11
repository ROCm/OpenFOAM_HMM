/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "alphaSgsJayatillekeWallFunctionFvPatchScalarField.H"
#include "LESModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar alphaSgsJayatillekeWallFunctionFvPatchScalarField::maxExp_ = 50.0;
scalar alphaSgsJayatillekeWallFunctionFvPatchScalarField::tolerance_ = 0.01;
label alphaSgsJayatillekeWallFunctionFvPatchScalarField::maxIters_ = 10;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void alphaSgsJayatillekeWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorIn
        (
            "alphaSgsJayatillekeWallFunctionFvPatchScalarField::checkType()"
        )
            << "Patch type for patch " << patch().name() << " must be wall\n"
            << "Current patch type is " << patch().type() << nl
            << exit(FatalError);
    }
}


scalar alphaSgsJayatillekeWallFunctionFvPatchScalarField::Psmooth
(
    const scalar Prat
) const
{
    return 9.24*(pow(Prat, 0.75) - 1.0)*(1.0 + 0.28*exp(-0.007*Prat));
}


scalar alphaSgsJayatillekeWallFunctionFvPatchScalarField::yPlusTherm
(
    const scalar P,
    const scalar Prat,
    const scalar E,
    const scalar kappa
) const
{
    scalar ypt = 11.0;

    for (int i=0; i<maxIters_; i++)
    {
        scalar f = ypt - (log(E*ypt)/kappa + P)/Prat;
        scalar df = 1.0 - 1.0/(ypt*kappa*Prat);
        scalar yptNew = ypt - f/df;

        if (yptNew < VSMALL)
        {
            return 0;
        }
        else if (mag(yptNew - ypt) < tolerance_)
        {
            return yptNew;
        }
        else
        {
            ypt = yptNew;
        }
     }

    return ypt;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphaSgsJayatillekeWallFunctionFvPatchScalarField::
alphaSgsJayatillekeWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{
    checkType();
}


alphaSgsJayatillekeWallFunctionFvPatchScalarField::
alphaSgsJayatillekeWallFunctionFvPatchScalarField
(
    const alphaSgsJayatillekeWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


alphaSgsJayatillekeWallFunctionFvPatchScalarField::
alphaSgsJayatillekeWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict)
{
    checkType();
}


alphaSgsJayatillekeWallFunctionFvPatchScalarField::
alphaSgsJayatillekeWallFunctionFvPatchScalarField
(
    const alphaSgsJayatillekeWallFunctionFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf)
{
    checkType();
}


alphaSgsJayatillekeWallFunctionFvPatchScalarField::
alphaSgsJayatillekeWallFunctionFvPatchScalarField
(
    const alphaSgsJayatillekeWallFunctionFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void alphaSgsJayatillekeWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    // Get info from the SGS model
    const LESModel& sgs = db().lookupObject<LESModel>("LESProperties");

    // Wall function constants
    const scalar E = sgs.E().value();
    const scalar kappa = sgs.kappa().value();
    const scalar Prt = sgs.Prt().value();

    // Field data
    const label patchI = patch().index();

    const scalarField& muw = sgs.mu().boundaryField()[patchI];
    const scalarField& muSgsw = sgs.muSgs()().boundaryField()[patchI];

    const scalarField& alphaw = sgs.alpha().boundaryField()[patchI];
    scalarField& alphaSgsw = *this;

    const fvPatchVectorField& Uw = sgs.U().boundaryField()[patchI];
    const scalarField magUp = mag(Uw.patchInternalField() - Uw);
    const scalarField magGradUw = mag(Uw.snGrad());

    const scalarField& rhow = sgs.rho().boundaryField()[patchI];
    const fvPatchScalarField& hw =
        patch().lookupPatchField<volScalarField, scalar>("h");

    const scalarField& ry = patch().deltaCoeffs();

    // Heat flux [W/m2] - lagging alphaSgsw
    const scalarField qDot = (alphaw + alphaSgsw)*hw.snGrad();

    // Populate boundary values
    forAll(alphaSgsw, faceI)
    {
        // Calculate uTau using Newton-Raphson iteration
        scalar uTau =
            sqrt((muSgsw[faceI] + muw[faceI])/rhow[faceI]*magGradUw[faceI]);

        if (uTau > ROOTVSMALL)
        {
            label iter = 0;
            scalar err = GREAT;

            do
            {
                scalar kUu = min(kappa*magUp[faceI]/uTau, maxExp_);
                scalar fkUu = exp(kUu) - 1.0 - kUu*(1.0 + 0.5*kUu);

                scalar f =
                    - uTau/(ry[faceI]*muw[faceI]/rhow[faceI])
                    + magUp[faceI]/uTau
                    + 1.0/E*(fkUu - 1.0/6.0*kUu*sqr(kUu));

                scalar df =
                    - 1.0/(ry[faceI]*muw[faceI]/rhow[faceI])
                    - magUp[faceI]/sqr(uTau)
                    - 1.0/E*kUu*fkUu/uTau;

                scalar uTauNew = uTau - f/df;
                err = mag((uTau - uTauNew)/uTau);
                uTau = uTauNew;

            } while (uTau>VSMALL && err>tolerance_ && ++iter<maxIters_);

            scalar yPlus = uTau/ry[faceI]/(muw[faceI]/rhow[faceI]);

            // Molecular Prandtl number
            scalar Pr = muw[faceI]/alphaw[faceI];

            // Molecular-to-turbulenbt Prandtl number ratio
            scalar Prat = Pr/Prt;

            // Thermal sublayer thickness
            scalar P = Psmooth(Prat);
            scalar yPlusTherm = this->yPlusTherm(P, Prat, E, kappa);

            // Evaluate new effective thermal diffusivity
            scalar alphaEff = 0.0;
            if (yPlus < yPlusTherm)
            {
                scalar A = qDot[faceI]*rhow[faceI]*uTau/ry[faceI];
                scalar B = qDot[faceI]*Pr*yPlus;
                scalar C = Pr*0.5*rhow[faceI]*uTau*sqr(magUp[faceI]);
                alphaEff = A/(B + C + VSMALL);
            }
            else
            {
                scalar A = qDot[faceI]*rhow[faceI]*uTau/ry[faceI];
                scalar B = qDot[faceI]*Prt*(1.0/kappa*log(E*yPlus) + P);
                scalar magUc = uTau/kappa*log(E*yPlusTherm) - mag(Uw[faceI]);
                scalar C =
                    0.5*rhow[faceI]*uTau
                   *(Prt*sqr(magUp[faceI]) + (Pr - Prt)*sqr(magUc));
                alphaEff = A/(B + C + VSMALL);
            }

            // Update turbulent thermal diffusivity
            alphaSgsw[faceI] = max(0.0, alphaEff - alphaw[faceI]);

            if (debug)
            {
                Info<< "    uTau           = " << uTau << nl
                    << "    Pr             = " << Pr << nl
                    << "    Prt            = " << Prt << nl
                    << "    qDot           = " << qDot[faceI] << nl
                    << "    yPlus          = " << yPlus << nl
                    << "    yPlusTherm     = " << yPlusTherm << nl
                    << "    alphaEff       = " << alphaEff << nl
                    << "    alphaw         = " << alphaw[faceI] << nl
                    << "    alphaSgsw      = " << alphaSgsw[faceI] << nl
                    << endl;
            }
        }
        else
        {
            alphaSgsw[faceI] = 0.0;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    alphaSgsJayatillekeWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
