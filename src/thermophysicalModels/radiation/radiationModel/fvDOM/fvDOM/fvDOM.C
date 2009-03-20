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

#include "fvDOM.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"

#include "absorptionEmissionModel.H"
#include "scatterModel.H"
#include "mathematicalConstants.H"
#include "radiationConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(fvDOM, 0);

        addToRunTimeSelectionTable
        (
            radiationModel,
            fvDOM,
            dictionary
        );
    }
}

//  Radiation solver iterator counter
Foam::label Foam::radiation::fvDOM::iterRadId = pTraits<label>::one;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::radiation::fvDOM::fvDOM(const volScalarField& T)
:
    radiationModel(typeName, T),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("G", dimMass/pow3(dimTime), 0.0)
    ),
    Qr_
    (
        IOobject
        (
            "Qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qr", dimMass/pow3(dimTime), 0.0)
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
    aj_(0),
    e_
    (
        IOobject
        (
            "e",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
    ),
    Ntheta_(readLabel(radiationModelCoeffs_.lookup("Ntheta"))),
    Nphi_(readLabel(radiationModelCoeffs_.lookup("Nphi"))),
    Ni_(0),
    lambdaj_(absorptionEmission_->nBands()),
    blackBody_(fileName("blackBodyEmissivePower"), "constant", lambdaj_, T)
{
    aj_.setSize(lambdaj_);
    if (mesh_.nSolutionD() == 3)    //3D
    {
        RadIntRay_.setSize(4.* Nphi_* Ntheta_);
        Ni_ = 4. * Nphi_ * Ntheta_;
        scalar deltaPhi = mathematicalConstant::pi / (2. * Nphi_);
        scalar deltaTheta = mathematicalConstant::pi / Ntheta_;
        label i = 0;
        for(label n = 1 ; n <= Ntheta_ ; n++)
        {
            for(label m = 1 ; m <= 4*Nphi_ ; m++)
            {
                scalar thetai = (2.*n - 1.)*deltaTheta/2.;
                scalar phii = (2.*m - 1.)*deltaPhi/2.;
                RadIntRay_.set
                (
                    i,
                    new radiativeIntensityRay
                    (
                        phii, thetai, deltaPhi,
                        deltaTheta, lambdaj_, mesh_,
                        absorptionEmission_, blackBody_
                    )
                );
                i++;
            }
        }
    }
    else
    {
        if (mesh_.nSolutionD() == 2)    //2D (X & Y)
        {
            scalar thetai = mathematicalConstant::pi/2.;
            scalar deltaTheta = mathematicalConstant::pi;
            RadIntRay_.setSize(4.* Nphi_);
            Ni_ = 4. * Nphi_;
            scalar deltaPhi = mathematicalConstant::pi / (2. * Nphi_);
            label i = 0;
            for(label m = 1 ; m <= 4*Nphi_ ; m++)
            {
                scalar phii = (2.*m - 1.)*deltaPhi/2.;
                RadIntRay_.set
                (
                    i,
                    new radiativeIntensityRay
                    (
                        phii, thetai, deltaPhi,
                        deltaTheta, lambdaj_, mesh_,
                        absorptionEmission_, blackBody_
                    )
                );
                i++;
            }
        }
        else    //1D (X)
        {
            scalar thetai = mathematicalConstant::pi/2.;
            scalar deltaTheta = mathematicalConstant::pi;
            RadIntRay_.setSize(2);
            Ni_ = 2.;
            scalar deltaPhi = mathematicalConstant::pi;
            label i = 0;
            for(label m = 1 ; m <= 2 ; m++)
            {
                scalar phii = (2*m - 1.)*deltaPhi/2.;
                RadIntRay_.set
                (
                    i,
                    new radiativeIntensityRay
                    (
                        phii, thetai, deltaPhi,
                        deltaTheta, lambdaj_, mesh_,
                        absorptionEmission_, blackBody_
                    )
                );
                i++;
            }

        }
    }


// Construct absorption for wave length
    for(label i=0; i < lambdaj_; i++)
    {
        volScalarField* volPtr= new volScalarField
        (
            IOobject
            (
                "aj_" + Foam::name(i) ,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            a_
        );

        aj_.set(i, volPtr);
    }

}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::fvDOM::~fvDOM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::fvDOM::read()
{
    if (radiationModel::read())
    {
        // nothing to read

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::radiation::fvDOM::correct()
{
    if (!radiation_ || !(iterRadId == nFlowIterPerRadIter_))
    {
        iterRadId++;
        return;
    }

    absorptionEmission_->correct(a_, aj_);

    updateBlackBodyEmission();

    scalar maxResidual = 0;
    scalar convergenceCriterion = 0;
    radiationModelCoeffs_.readIfPresent("convergence", convergenceCriterion);
    label radIter = 0;
    do
    {
         radIter ++;
        for(label i = 0; i < Ni_; i++) //
        {
            maxResidual = 0;
            scalar maxBandResidual = RadIntRay_[i].correct(this);
            maxResidual = max(maxBandResidual, maxResidual);
        }

        Info << "Radiation solver Iter: " <<  radIter << endl;

    }while(maxResidual > convergenceCriterion);

    updateG();

    iterRadId = pTraits<label>::one;
}

Foam::tmp<Foam::volScalarField> Foam::radiation::fvDOM::Rp() const
{

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Rp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            4.0*a_*radiation::sigmaSB //absorptionEmission_->a()
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::radiation::fvDOM::Ru() const
{

    const DimensionedField<scalar, volMesh>& G =
        G_.dimensionedInternalField();
    const DimensionedField<scalar, volMesh> E =
        absorptionEmission_->ECont()().dimensionedInternalField();
    const DimensionedField<scalar, volMesh> a =
        a_.dimensionedInternalField(); //absorptionEmission_->aCont()()

    return  a*G - 4.0*E;
}

void Foam::radiation::fvDOM::updateBlackBodyEmission()
{

    for(label j=0; j < lambdaj_; j++)
    {
          blackBody_.correct(j, absorptionEmission_->bands(j));
    }

}

void Foam::radiation::fvDOM::updateG()
{

    G_ = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
    Qr_ = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);

    for(label i = 0; i < Ni_; i++)
    {
        RadIntRay_[i].addIntensity();
        G_ +=  RadIntRay_[i].I()* RadIntRay_[i].omegai();
        Qr_ += RadIntRay_[i].Qri();
    }
}
// ************************************************************************* //
