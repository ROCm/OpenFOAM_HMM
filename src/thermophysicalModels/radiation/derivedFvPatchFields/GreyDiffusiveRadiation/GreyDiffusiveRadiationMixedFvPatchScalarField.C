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

#include "GreyDiffusiveRadiationMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

#include "radiationModel.H"
#include "fvDOM.H"
#include "radiationConstants.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::GreyDiffusiveRadiationMixedFvPatchField::
GreyDiffusiveRadiationMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    TName_("undefined"),
    emissivity_(0.0),
    myRayIndex_(0),
    myWaveLengthIndex_(0),
    myRayIsInit_(-1),
    qr_(0)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
    qr_.setSize(p.size());
}


Foam::radiation::GreyDiffusiveRadiationMixedFvPatchField::
GreyDiffusiveRadiationMixedFvPatchField
(
    const GreyDiffusiveRadiationMixedFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    emissivity_(ptf.emissivity_),
    myRayIndex_(ptf.myRayIndex_),
    myWaveLengthIndex_(ptf.myWaveLengthIndex_),
    myRayIsInit_(ptf.myRayIsInit_),
    qr_(ptf.qr_)
{}


Foam::radiation::GreyDiffusiveRadiationMixedFvPatchField::
GreyDiffusiveRadiationMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    TName_(dict.lookup("T")),
    emissivity_(readScalar(dict.lookup("emissivity"))),
    myRayIndex_(0),
    myWaveLengthIndex_(0),
    myRayIsInit_(-1),
    qr_(0)
{
    const scalarField& Tp =
        patch().lookupPatchField<volScalarField, scalar>(TName_);

    refValue() = emissivity_*4.0*radiation::sigmaSB.value()*pow4(Tp) /
                 Foam::mathematicalConstant::pi;
    refGrad() = 0.0;

    qr_.setSize(p.size());

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


Foam::radiation::GreyDiffusiveRadiationMixedFvPatchField::
GreyDiffusiveRadiationMixedFvPatchField
(
    const GreyDiffusiveRadiationMixedFvPatchField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    emissivity_(ptf.emissivity_),
    myRayIndex_(ptf.myRayIndex_),
    myWaveLengthIndex_(ptf.myWaveLengthIndex_),
    myRayIsInit_(ptf.myRayIsInit_),
    qr_(ptf.qr_)
{}


Foam::radiation::GreyDiffusiveRadiationMixedFvPatchField::
GreyDiffusiveRadiationMixedFvPatchField
(
    const GreyDiffusiveRadiationMixedFvPatchField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    emissivity_(ptf.emissivity_),
    myRayIndex_(ptf.myRayIndex_),
    myWaveLengthIndex_(ptf.myWaveLengthIndex_),
    myRayIsInit_(ptf.myRayIsInit_),
    qr_(ptf.qr_)

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::GreyDiffusiveRadiationMixedFvPatchField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    scalarField::autoMap(m);

}


void Foam::radiation::GreyDiffusiveRadiationMixedFvPatchField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

//    const GreyDiffusiveRadiationMixedFvPatchField& mrptf =
        refCast<const GreyDiffusiveRadiationMixedFvPatchField>(ptf);
}


void Foam::radiation::GreyDiffusiveRadiationMixedFvPatchField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalarField& Tp =
        patch().lookupPatchField<volScalarField, scalar>(TName_);

    const radiationModel& rad =
            db().lookupObject<radiationModel>("radiationProperties");

    const fvDOM& Dom(refCast<const fvDOM>(rad));

    const label patchi = patch().index();

    if(Dom.lambdaj() == 1)
    {
        if (myRayIsInit_ == -1)
        {
            for(label i=0; i < Dom.Ni() ; i++)
            {
                 for(label j=0; j < Dom.lambdaj() ; j++)
                {
                    const volScalarField& radiationField =
                                        Dom.RadIntRayiLambdaj(i,j);
                    if (&(radiationField.internalField()) ==
                        &dimensionedInternalField())
                    {
                        myRayIndex_ = i;
                        myWaveLengthIndex_ = j;
                        myRayIsInit_ = 0.;
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
            "GreyDiffusiveRadiationMixedFvPatchField::"
            "updateCoeffs"
        )   << " a grey boundary condition is used with a non-grey"
            << "absorption model"
            << exit(FatalError);
    }

    vectorField n = patch().Sf()/patch().magSf();

    scalarField& Iw = *(this);

    qr_ =  Iw *(-n & Dom.RadIntRay(myRayIndex_).Di());

    Dom.RadIntRay(myRayIndex_).add(qr_,patchi);

    forAll(Iw, faceI)
    {

        scalar Ir = 0.0;

        for(label i=0; i < Dom.Ni() ; i++) //
        {
            const vector& si = Dom.RadIntRay(i).Si();

            const scalarField&  Iface = Dom.RadIntRay(i).Ilambdaj
            (
                myWaveLengthIndex_
            ).boundaryField()[patch().index()];

            scalar InOut = -n[faceI] & si;

            if (InOut < 0.) // qin into the wall
            {
                const vector& di = Dom.RadIntRay(i).Di();
                Ir += Iface[faceI]*mag(n[faceI] & di);
            }
        }

        const vector& mySi = Dom.RadIntRay(myRayIndex_).Si();

        scalar InOut = -n[faceI] & mySi;

        if (InOut > 0.) //direction out of the wall
        {
            refGrad()[faceI] = 0.0;
            valueFraction()[faceI] = 1.0;
            refValue()[faceI] = ((1. - emissivity_) * Ir +
                    emissivity_*radiation::sigmaSB.value()*pow4(Tp[faceI])) /
                    Foam::mathematicalConstant::pi;

        }
        else if (InOut < 0.) //direction into the wall
        {
            valueFraction()[faceI] = 0.0;
            refGrad()[faceI] = 0.0;
            refValue()[faceI] = 0.0; //not used
        }
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::radiation::GreyDiffusiveRadiationMixedFvPatchField::write(Ostream&
os) const
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
        GreyDiffusiveRadiationMixedFvPatchField
    );
}
}


// ************************************************************************* //
