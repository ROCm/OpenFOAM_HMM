/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2022 Upstream CFD GmbH
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "SpalartAllmarasDDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDDES<BasicTurbulenceModel>::fd
(
    const volScalarField& magGradU
) const
{
    return
        1
      - tanh(pow(this->Cd1_*this->r(this->nuEff(), magGradU, this->y_), Cd2_));
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDDES<BasicTurbulenceModel>::Stilda
(
    const volScalarField& chi,
    const volScalarField& fv1,
    const volTensorField& gradU,
    const volScalarField& dTilda
) const
{
    if (this->useSigma_)
    {
        const volScalarField& lRAS(this->y_);
        const volScalarField fv2(this->fv2(chi, fv1));
        const volScalarField lLES(this->lengthScaleLES(chi, fv1));
        const volScalarField Omega(this->Omega(gradU));
        const volScalarField Ssigma(this->Ssigma(gradU));
        const volScalarField SsigmaDES
        (
            Omega - fd(mag(gradU))*pos(lRAS - lLES)*(Omega - Ssigma)
        );

        return
            max
            (
                SsigmaDES + fv2*this->nuTilda_/sqr(this->kappa_*dTilda),
                this->Cs_*SsigmaDES
            );
    }

    return
        SpalartAllmarasBase<DESModel<BasicTurbulenceModel>>::Stilda
        (
            chi,
            fv1,
            gradU,
            dTilda
        );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDDES<BasicTurbulenceModel>::dTilda
(
    const volScalarField& chi,
    const volScalarField& fv1,
    const volTensorField& gradU
) const
{
    const volScalarField& lRAS(this->y_);
    const volScalarField lLES(this->lengthScaleLES(chi, fv1));
    const dimensionedScalar l0(dimLength, Zero);

    return max
    (
        lRAS - fd(mag(gradU))*max(lRAS - lLES, l0),
        dimensionedScalar("small", dimLength, SMALL)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
SpalartAllmarasDDES<BasicTurbulenceModel>::SpalartAllmarasDDES
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
    SpalartAllmarasDES<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    ),

    Cd1_
    (
        this->useSigma_
          ? dimensioned<scalar>::getOrAddToDict
            (
                "Cd1Sigma",
                this->coeffDict_,
                10
            )
          : dimensioned<scalar>::getOrAddToDict
            (
                "Cd1",
                this->coeffDict_,
                8
            )
    ),
    Cd2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cd2",
            this->coeffDict_,
            3
        )
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool SpalartAllmarasDDES<BasicTurbulenceModel>::read()
{
    if (SpalartAllmarasDES<BasicTurbulenceModel>::read())
    {
        Cd1_.readIfPresent(this->coeffDict());
        Cd2_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDDES<BasicTurbulenceModel>::fd() const
{
    return fd(mag(fvc::grad(this->U_)));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
