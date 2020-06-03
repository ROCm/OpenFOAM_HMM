/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "SpalartAllmarasDES.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::chi() const
{
    return nuTilda_/this->nu();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::fv1
(
    const volScalarField& chi
) const
{
    const volScalarField chi3("chi3", pow3(chi));
    return chi3/(chi3 + pow3(Cv1_));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::fv2
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    return 1.0 - chi/(1.0 + chi*fv1);
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::ft2
(
    const volScalarField& chi
) const
{
    return Ct3_*exp(-Ct4_*sqr(chi));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::Omega
(
    const volTensorField& gradU
) const
{
    return sqrt(2.0)*mag(skew(gradU));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::Stilda
(
    const volScalarField& chi,
    const volScalarField& fv1,
    const volScalarField& Omega,
    const volScalarField& dTilda
) const
{
    return
    (
        max
        (
            Omega
          + fv2(chi, fv1)*nuTilda_/sqr(kappa_*dTilda),
            Cs_*Omega
        )
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::r
(
    const volScalarField& nur,
    const volScalarField& Omega,
    const volScalarField& dTilda
) const
{
    tmp<volScalarField> tr
    (
        min
        (
            nur
           /(
                max
                (
                    Omega,
                    dimensionedScalar("SMALL", Omega.dimensions(), SMALL)
                )
                *sqr(kappa_*dTilda)
            ),
            scalar(10)
        )
    );
    tr.ref().boundaryFieldRef() == 0.0;

    return tr;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::fw
(
    const volScalarField& Omega,
    const volScalarField& dTilda
) const
{
    const volScalarField r(this->r(nuTilda_, Omega, dTilda));
    const volScalarField g(r + Cw2_*(pow6(r) - r));

    return g*pow((1 + pow6(Cw3_))/(pow6(g) + pow6(Cw3_)), 1.0/6.0);
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::psi
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    tmp<volScalarField> tpsi
    (
        new volScalarField
        (
            IOobject
            (
                type() + ":psi",
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("one", dimless, 1)
        )
    );

    if (lowReCorrection_)
    {
        volScalarField& psi = tpsi.ref();

        const volScalarField fv2(this->fv2(chi, fv1));
        const volScalarField ft2(this->ft2(chi));

        psi =
            sqrt
            (
                min
                (
                    scalar(100),
                    (1 - Cb1_/(Cw1_*sqr(kappa_)*fwStar_)*(ft2 + (1 - ft2)*fv2))
                   /max(SMALL, (fv1*max(scalar(1e-10), 1 - ft2)))
                )
            );
    }

    return tpsi;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::dTilda
(
    const volScalarField& chi,
    const volScalarField& fv1,
    const volTensorField& gradU
) const
{
    tmp<volScalarField> tdTilda(psi(chi, fv1)*CDES_*this->delta());
    min(tdTilda.ref().ref(), tdTilda(), y_);
    return tdTilda;
}


template<class BasicTurbulenceModel>
void SpalartAllmarasDES<BasicTurbulenceModel>::correctNut
(
    const volScalarField& fv1
)
{
    this->nut_ = nuTilda_*fv1;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void SpalartAllmarasDES<BasicTurbulenceModel>::correctNut()
{
    correctNut(fv1(this->chi()));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
SpalartAllmarasDES<BasicTurbulenceModel>::SpalartAllmarasDES
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    DESModel<BasicTurbulenceModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    sigmaNut_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaNut",
            this->coeffDict_,
            0.66666
        )
    ),
    kappa_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    Cb1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cb1",
            this->coeffDict_,
            0.1355
        )
    ),
    Cb2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cb2",
            this->coeffDict_,
            0.622
        )
    ),
    Cw1_(Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_),
    Cw2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cw2",
            this->coeffDict_,
            0.3
        )
    ),
    Cw3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cw3",
            this->coeffDict_,
            2.0
        )
    ),
    Cv1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cv1",
            this->coeffDict_,
            7.1
        )
    ),
    Cs_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            0.3
        )
    ),
    CDES_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "CDES",
            this->coeffDict_,
            0.65
        )
    ),
    ck_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "ck",
            this->coeffDict_,
            0.07
        )
    ),
    lowReCorrection_
    (
        Switch::getOrAddToDict
        (
            "lowReCorrection",
            this->coeffDict_,
            true
        )
    ),
    Ct3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ct3",
            this->coeffDict_,
            1.2
        )
    ),
    Ct4_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ct4",
            this->coeffDict_,
            0.5
        )
    ),
    fwStar_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "fwStar",
            this->coeffDict_,
            0.424
        )
    ),

    nuTilda_
    (
        IOobject
        (
            "nuTilda",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    y_(wallDist::New(this->mesh_).y())
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool SpalartAllmarasDES<BasicTurbulenceModel>::read()
{
    if (DESModel<BasicTurbulenceModel>::read())
    {
        sigmaNut_.readIfPresent(this->coeffDict());
        kappa_.readIfPresent(*this);

        Cb1_.readIfPresent(this->coeffDict());
        Cb2_.readIfPresent(this->coeffDict());
        Cw1_ = Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_;
        Cw2_.readIfPresent(this->coeffDict());
        Cw3_.readIfPresent(this->coeffDict());
        Cv1_.readIfPresent(this->coeffDict());
        Cs_.readIfPresent(this->coeffDict());

        CDES_.readIfPresent(this->coeffDict());
        ck_.readIfPresent(this->coeffDict());

        lowReCorrection_.readIfPresent("lowReCorrection", this->coeffDict());
        Ct3_.readIfPresent(this->coeffDict());
        Ct4_.readIfPresent(this->coeffDict());
        fwStar_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::
DnuTildaEff() const
{
    return tmp<volScalarField>
    (
        new volScalarField("DnuTildaEff", (nuTilda_ + this->nu())/sigmaNut_)
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::k() const
{
    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));

    tmp<volScalarField> tdTilda
    (
        new volScalarField
        (
            IOobject
            (
                "dTilda",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            dimensionedScalar(dimLength, Zero),
            zeroGradientFvPatchField<vector>::typeName
        )
    );
    volScalarField& dTildaF = tdTilda.ref();
    dTildaF = dTilda(chi, fv1, fvc::grad(this->U_));
    dTildaF.correctBoundaryConditions();

    return sqr(this->nut()/ck_/dTildaF);
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::LESRegion() const
{
    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));

    tmp<volScalarField> tLESRegion
    (
        new volScalarField
        (
            IOobject
            (
                "DES::LESRegion",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            neg(dTilda(chi, fv1, fvc::grad(this->U_)) - y_)
        )
    );

    return tLESRegion;
}


template<class BasicTurbulenceModel>
void SpalartAllmarasDES<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    DESModel<BasicTurbulenceModel>::correct();

    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));
    const volScalarField ft2(this->ft2(chi));

    tmp<volTensorField> tgradU = fvc::grad(U);
    const volScalarField Omega(this->Omega(tgradU()));
    const volScalarField dTilda(this->dTilda(chi, fv1, tgradU()));
    const volScalarField Stilda(this->Stilda(chi, fv1, Omega, dTilda));

    tmp<fvScalarMatrix> nuTildaEqn
    (
        fvm::ddt(alpha, rho, nuTilda_)
      + fvm::div(alphaRhoPhi, nuTilda_)
      - fvm::laplacian(alpha*rho*DnuTildaEff(), nuTilda_)
      - Cb2_/sigmaNut_*alpha*rho*magSqr(fvc::grad(nuTilda_))
     ==
        Cb1_*alpha*rho*Stilda*nuTilda_*(scalar(1) - ft2)
      - fvm::Sp
        (
            (Cw1_*fw(Stilda, dTilda) - Cb1_/sqr(kappa_)*ft2)
           *alpha*rho*nuTilda_/sqr(dTilda),
            nuTilda_
        )
      + fvOptions(alpha, rho, nuTilda_)
    );

    nuTildaEqn.ref().relax();
    fvOptions.constrain(nuTildaEqn.ref());
    solve(nuTildaEqn);
    fvOptions.correct(nuTilda_);
    bound(nuTilda_, dimensionedScalar(nuTilda_.dimensions(), Zero));
    nuTilda_.correctBoundaryConditions();

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
