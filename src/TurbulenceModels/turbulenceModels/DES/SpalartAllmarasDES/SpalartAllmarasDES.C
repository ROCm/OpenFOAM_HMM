/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::psi
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    auto tpsi = tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::scopedName(type(), "psi"),
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless, Zero)
    );

    if (lowReCorrection_)
    {
        auto& psi = tpsi.ref();

        const volScalarField fv2(this->fv2(chi, fv1));
        const volScalarField ft2(this->ft2(chi));

        psi =
            sqrt
            (
                min
                (
                    scalar(100),
                    (1 - this->Cb1_/(this->Cw1_*sqr(this->kappa_)*fwStar_)
                   *(ft2 + (1 - ft2)*fv2))
                   /max(SMALL, (fv1*max(scalar(1e-10), 1 - ft2)))
                )
            );
    }

    return tpsi;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::lengthScaleLES
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    return psi(chi, fv1)*CDES_*this->delta();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::dTilda
(
    const volScalarField& chi,
    const volScalarField& fv1,
    const volTensorField& gradU
) const
{
    // Initialise with LES length scale
    tmp<volScalarField> tdTilda = lengthScaleLES(chi, fv1);

    // Take min vs wall distance
    min(tdTilda.ref(), tdTilda(), this->y_);
    return tdTilda;
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
    SpalartAllmarasBase<DESModel<BasicTurbulenceModel>>
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

    CDES_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "CDES",
            this->coeffDict_,
            0.65
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
    fwStar_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "fwStar",
            this->coeffDict_,
            0.424
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
bool SpalartAllmarasDES<BasicTurbulenceModel>::read()
{
    if (SpalartAllmarasBase<DESModel<BasicTurbulenceModel>>::read())
    {
        CDES_.readIfPresent(this->coeffDict());
        lowReCorrection_.readIfPresent("lowReCorrection", this->coeffDict());
        fwStar_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
tmp<volScalarField>
SpalartAllmarasDES<BasicTurbulenceModel>::LESRegion() const
{
    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));

    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::scopedName("DES", "LESRegion"),
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        neg(dTilda(chi, fv1, fvc::grad(this->U_)) - this->y_)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
