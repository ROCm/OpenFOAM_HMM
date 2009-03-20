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
    nTheta_(readLabel(coeffs_.lookup("nTheta"))),
    nPhi_(readLabel(coeffs_.lookup("nPhi"))),
    nRay_(0),
    nLambda_(absorptionEmission_->nBands()),
    aj_(nLambda_),
    blackBody_
    (
        fileName("blackBodyEmissivePower"),
        mesh_.time().constant(),
        nLambda_,
        T
    ),
    IRay_(0),
    convergence_(coeffs_.lookupOrDefault<scalar>("convergence", 0.0))
{
    if (mesh_.nSolutionD() == 3)    //3D
    {
        IRay_.setSize(4*nPhi_*nTheta_);
        nRay_ = 4.0*nPhi_*nTheta_;
        scalar deltaPhi = mathematicalConstant::pi/(2.0*nPhi_);
        scalar deltaTheta = mathematicalConstant::pi/nTheta_;
        label i = 0;
        for (label n = 1; n <= nTheta_; n++)
        {
            for (label m = 1; m <= 4*nPhi_; m++)
            {
                scalar thetai = (2.0*n - 1.0)*deltaTheta/2.0;
                scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
                IRay_.set
                (
                    i,
                    new radiativeIntensityRay
                    (
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        nLambda_,
                        mesh_,
                        absorptionEmission_,
                        blackBody_
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
            scalar thetai = mathematicalConstant::pi/2.0;
            scalar deltaTheta = mathematicalConstant::pi;
            IRay_.setSize(4*nPhi_);
            nRay_ = 4.0*nPhi_;
            scalar deltaPhi = mathematicalConstant::pi/(2.0*nPhi_);
            label i = 0;
            for (label m = 1; m <= 4*nPhi_; m++)
            {
                scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
                IRay_.set
                (
                    i,
                    new radiativeIntensityRay
                    (
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        nLambda_,
                        mesh_,
                        absorptionEmission_,
                        blackBody_
                    )
                );
                i++;
            }
        }
        else    //1D (X)
        {
            scalar thetai = mathematicalConstant::pi/2.0;
            scalar deltaTheta = mathematicalConstant::pi;
            IRay_.setSize(2);
            nRay_ = 2.0;
            scalar deltaPhi = mathematicalConstant::pi;
            label i = 0;
            for (label m = 1; m <= 2; m++)
            {
                scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
                IRay_.set
                (
                    i,
                    new radiativeIntensityRay
                    (
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        nLambda_,
                        mesh_,
                        absorptionEmission_,
                        blackBody_
                    )
                );
                i++;
            }

        }
    }


    // Construct absorption field for each wavelength
    forAll(aj_, lambdaI)
    {
        aj_.set
        (
            lambdaI,
            new volScalarField
            (
                IOobject
                (
                    "aj_" + Foam::name(lambdaI) ,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                a_
            )
        );
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

//        coeffs_.lookup("nTheta") >> nTheta_;
//        coeffs_.lookup("nPhi") >> nPhi_;

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
    label radIter = 0;
    do
    {
        radIter ++;
        forAll(IRay_, rayI)
        {
            maxResidual = 0.0;
            scalar maxBandResidual = IRay_[rayI].correct(this);
            maxResidual = max(maxBandResidual, maxResidual);
        }

        Info << "Radiation solver Iter: " <<  radIter << endl;

    } while(maxResidual > convergence_);

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
    for (label j=0; j < nLambda_; j++)
    {
        blackBody_.correct(j, absorptionEmission_->bands(j));
    }
}


void Foam::radiation::fvDOM::updateG()
{
    G_ = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
    Qr_ = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);

    forAll(IRay_, rayI)
    {
        IRay_[rayI].addIntensity();
        G_ += IRay_[rayI].I()*IRay_[rayI].omega();
        Qr_ += IRay_[rayI].Qr();
    }
}


// ************************************************************************* //
