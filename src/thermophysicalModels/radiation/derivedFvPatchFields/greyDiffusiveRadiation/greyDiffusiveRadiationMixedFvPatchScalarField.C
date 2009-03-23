/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "greyDiffusiveRadiationMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

#include "fvDOM.H"
#include "radiationConstants.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::
greyDiffusiveRadiationMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    TName_("undefinedT"),
    emissivity_(0.0),
    rayId_(0),
    lambdaId_(0)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::
greyDiffusiveRadiationMixedFvPatchScalarField
(
    const greyDiffusiveRadiationMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    emissivity_(ptf.emissivity_),
    rayId_(ptf.rayId_),
    lambdaId_(ptf.lambdaId_)
{}


Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::
greyDiffusiveRadiationMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    TName_(dict.lookup("T")),
    emissivity_(readScalar(dict.lookup("emissivity"))),
    rayId_(-1),
    lambdaId_(-1)
{
    const scalarField& Tp =
        patch().lookupPatchField<volScalarField, scalar>(TName_);

    refValue() =
        emissivity_*4.0*radiation::sigmaSB.value()*pow4(Tp)
       /Foam::mathematicalConstant::pi;

    refGrad() = 0.0;

    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchScalarField::operator=(refValue());
    }
}


Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::
greyDiffusiveRadiationMixedFvPatchScalarField
(
    const greyDiffusiveRadiationMixedFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    emissivity_(ptf.emissivity_),
    rayId_(ptf.rayId_),
    lambdaId_(ptf.lambdaId_)
{}


Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::
greyDiffusiveRadiationMixedFvPatchScalarField
(
    const greyDiffusiveRadiationMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    emissivity_(ptf.emissivity_),
    rayId_(ptf.rayId_),
    lambdaId_(ptf.lambdaId_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::
updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalarField& Tp =
        patch().lookupPatchField<volScalarField, scalar>(TName_);

    const fvDOM& dom = db().lookupObject<fvDOM>("radiationProperties");

    const label patchI = patch().index();

    if (dom.nLambda() == 1)
    {
        if (rayId_ == -1)
        {
            for (label rayI=0; rayI < dom.nRay(); rayI++)
            {
                for (label lambdaI=0; lambdaI < dom.nLambda(); lambdaI++)
                {
                    const volScalarField& radiationField =
                        dom.IRayLambda(rayI, lambdaI);
                    if
                    (
                        &(radiationField.internalField())
                     == &dimensionedInternalField()
                    )
                    {
                        rayId_ = rayI;
                        lambdaId_ = lambdaI;
                        break;
                    }
                }
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "Foam::radiation::"
            "greyDiffusiveRadiationMixedFvPatchScalarField::"
            "updateCoeffs"
        )   << " a grey boundary condition is used with a non-grey"
            << "absorption model"
            << exit(FatalError);
    }

    scalarField& Iw = *this;
    vectorField n = patch().Sf()/patch().magSf();

    radiativeIntensityRay& ray =
        const_cast<radiativeIntensityRay&>(dom.IRay(rayId_));

    ray.Qr().boundaryField()[patchI] += Iw*(-n & ray.dAve());

    forAll(Iw, faceI)
    {
        scalar Ir = 0.0;

        for (label rayI=0; rayI < dom.nRay(); rayI++)
        {
            const vector& d = dom.IRay(rayI).d();

            const scalarField& Iface =
                dom.IRay(rayI).ILambda(lambdaId_).boundaryField()[patchI];

            if ((-n[faceI] & d) < 0.0) // qin into the wall
            {
                const vector& dAve = dom.IRay(rayI).dAve();
                Ir += Iface[faceI]*mag(n[faceI] & dAve);
            }
        }

        const vector& d = dom.IRay(rayId_).d();

        if ((-n[faceI] & d) > 0.) //direction out of the wall
        {
            refGrad()[faceI] = 0.0;
            valueFraction()[faceI] = 1.0;
            refValue()[faceI] =
                (
                    Ir*(1.0 - emissivity_)
                  + emissivity_*radiation::sigmaSB.value()*pow4(Tp[faceI])
                )
               /mathematicalConstant::pi;

        }
        else //direction into the wall
        {
            valueFraction()[faceI] = 0.0;
            refGrad()[faceI] = 0.0;
            refValue()[faceI] = 0.0; //not used
        }
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::radiation::greyDiffusiveRadiationMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("T") << TName_ << token::END_STATEMENT << nl;
    os.writeKeyword("emissivity") << emissivity_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
    makePatchTypeField
    (
        fvPatchScalarField,
        greyDiffusiveRadiationMixedFvPatchScalarField
    );
}
}


// ************************************************************************* //
