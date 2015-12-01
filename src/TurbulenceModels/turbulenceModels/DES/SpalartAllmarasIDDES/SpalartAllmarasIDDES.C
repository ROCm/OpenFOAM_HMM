/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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

#include "SpalartAllmarasIDDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasIDDES<BasicTurbulenceModel>::alpha() const
{
    // Equation 9 (plus limits)
    return max
    (
        0.25 - this->y_/static_cast<const volScalarField&>(IDDESDelta_.hmax()),
        scalar(-5)
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasIDDES<BasicTurbulenceModel>::ft
(
    const volScalarField& magGradU
) const
{
    // Equation 13
    return tanh(pow3(sqr(ct_)*rd(this->nut_, magGradU)));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasIDDES<BasicTurbulenceModel>::fl
(
    const volScalarField& magGradU
) const
{
    // Equation 13
    return tanh(pow(sqr(cl_)*rd(this->nu(), magGradU), 10));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasIDDES<BasicTurbulenceModel>::rd
(
    const volScalarField& nur,
    const volScalarField& magGradU
) const
{
    return min
    (
        nur
       /(
           max
           (
               magGradU,
               dimensionedScalar("SMALL", magGradU.dimensions(), SMALL)
           )*sqr(this->kappa_*this->y_)
       ),
       scalar(10)
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasIDDES<BasicTurbulenceModel>::fdt
(
    const volScalarField& magGradU
) const
{
    // Related to equation 16
    return 1 - tanh(pow3(8*rd(this->nuEff(), magGradU)));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasIDDES<BasicTurbulenceModel>::dTilda
(
    const volScalarField& chi,
    const volScalarField& fv1,
    const volTensorField& gradU
) const
{
    const volScalarField alpha(this->alpha());
    const volScalarField expTerm(exp(sqr(alpha)));
    const volScalarField magGradU(mag(gradU));

    // Equation 9
    tmp<volScalarField> fB = min(2*pow(expTerm, -9.0), scalar(1));

    // Equation 11
    tmp<volScalarField> fe1 =
        2*(pos(alpha)*pow(expTerm, -11.09) + neg(alpha)*pow(expTerm, -9.0));

    // Equation 12
    tmp<volScalarField> fe2 = 1 - max(ft(magGradU), fl(magGradU));

    // Equation 10
    const volScalarField psi(this->psi(chi, fv1));
    tmp<volScalarField> fe = max(fe1 - 1, scalar(0))*psi*fe2;

    // Equation 16
    const volScalarField fdTilda(max(1 - fdt(magGradU), fB));

    // Equation 17 (plus limits)
    return max
    (
        fdTilda*(1 + fe)*this->y_ + (1 - fdTilda)*psi*this->CDES_*this->delta(),
        dimensionedScalar("SMALL", dimLength, SMALL)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
SpalartAllmarasIDDES<BasicTurbulenceModel>::SpalartAllmarasIDDES
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
        propertiesName
    ),
    fwStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "fwStar",
            this->coeffDict_,
            0.424
        )
    ),
    cl_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cl",
            this->coeffDict_,
            3.55
        )
    ),
    ct_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ct",
            this->coeffDict_,
            1.63
        )
    ),
    IDDESDelta_(refCast<IDDESDelta>(this->delta_()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool SpalartAllmarasIDDES<BasicTurbulenceModel>::read()
{
    if (SpalartAllmarasDES<BasicTurbulenceModel>::read())
    {
        fwStar_.readIfPresent(this->coeffDict());
        cl_.readIfPresent(this->coeffDict());
        ct_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
