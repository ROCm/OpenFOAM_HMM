/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
tmp<volScalarField> SpalartAllmarasDDES<BasicTurbulenceModel>::rd
(
    const volScalarField& magGradU
) const
{
    tmp<volScalarField> tr
    (
        min
        (
            this->nuEff()
           /(
                max
                (
                    magGradU,
                    dimensionedScalar("SMALL", magGradU.dimensions(), SMALL)
                )
               *sqr(this->kappa_*this->y_)
            ),
            scalar(10)
        )
    );
    tr.ref().boundaryFieldRef() == 0.0;

    return tr;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDDES<BasicTurbulenceModel>::fd
(
    const volScalarField& magGradU
) const
{
    return 1 - tanh(pow(Cd1_*rd(magGradU), Cd2_));
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDDES<BasicTurbulenceModel>::dTilda
(
    const volScalarField& chi,
    const volScalarField& fv1,
    const volTensorField& gradU
) const
{
    const volScalarField& lRAS(this->y_);
    const volScalarField lLES(this->psi(chi, fv1)*this->CDES_*this->delta());

    return max
    (
        lRAS
      - fd(mag(gradU))
       *max
        (
            lRAS - lLES,
            dimensionedScalar(dimLength, Zero)
        ),
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
        dimensioned<scalar>::getOrAddToDict
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
