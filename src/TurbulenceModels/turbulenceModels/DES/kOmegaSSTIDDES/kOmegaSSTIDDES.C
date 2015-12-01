/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "kOmegaSSTIDDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
const IDDESDelta& kOmegaSSTIDDES<BasicTurbulenceModel>::setDelta() const
{
    if (!isA<IDDESDelta>(this->delta_()))
    {
        FatalErrorIn
        (
            "const kOmegaSSTIDDES<BasicTurbulenceModel>::setDelta() const"
        )
            << "The delta function must be set to a " << IDDESDelta::typeName
            << " -based model" << exit(FatalError);
    }

    return refCast<const IDDESDelta>(this->delta_());
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTIDDES<BasicTurbulenceModel>::alpha() const
{
    return max(0.25 - this->y_/IDDESDelta_.hmax(), scalar(-5));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTIDDES<BasicTurbulenceModel>::ft
(
    const volScalarField& magGradU
) const
{
    return tanh(pow3(sqr(ct_)*rd(this->nut_, magGradU)));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTIDDES<BasicTurbulenceModel>::fl
(
    const volScalarField& magGradU
) const
{
    return tanh(pow(sqr(cl_)*rd(this->nu(), magGradU), 10));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTIDDES<BasicTurbulenceModel>::rd
(
    const volScalarField& nur,
    const volScalarField& magGradU
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
                    magGradU,
                    dimensionedScalar("SMALL", magGradU.dimensions(), SMALL)
                )
                *sqr(this->kappa_*this->y_)
            ),
            scalar(10)
        )
    );
    tr().boundaryField() == 0.0;

    return tr;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTIDDES<BasicTurbulenceModel>::fd
(
    const volScalarField& magGradU
) const
{
    return 1 - tanh(pow(cdt1_*rd(this->nuEff(), magGradU), cdt2_));
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTIDDES<BasicTurbulenceModel>::dTilda
(
    const volScalarField& magGradU,
    const volScalarField& CDES
) const
{
    const volScalarField& k = this->k_;
    const volScalarField& omega = this->omega_;

    const volScalarField lRAS(sqrt(k)/(this->betaStar_*omega));
    const volScalarField lLES(CDES*this->delta());
    const dimensionedScalar d0("SMALL", dimLength, SMALL);

    const volScalarField alpha(this->alpha());
    const volScalarField expTerm(exp(sqr(alpha)));

    tmp<volScalarField> fStep = min(2*pow(expTerm, -9.0), scalar(1));
    const volScalarField fHyb(max(1 - fd(magGradU), fStep));
    // Simplified version where fRestore = 0
    // return max(d0, fHyb*lRAS + (1 - fHyb)*lLES);

    // Original form
    tmp<volScalarField> fHill =
        2*(pos(alpha)*pow(expTerm, -11.09) + neg(alpha)*pow(expTerm, -9.0));
    tmp<volScalarField> fAmp = 1 - max(ft(magGradU), fl(magGradU));
    tmp<volScalarField> fRestore = max(fHill - 1, scalar(0))*fAmp;
    return max(d0, fHyb*(1 + fRestore)*lRAS + (1 - fHyb)*lLES);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaSSTIDDES<BasicTurbulenceModel>::kOmegaSSTIDDES
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
    kOmegaSSTDES<BasicTurbulenceModel>
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
    cdt1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cdt1",
            this->coeffDict_,
            20
        )
    ),
    cdt2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cdt2",
            this->coeffDict_,
            3
        )
    ),
    cl_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cl",
            this->coeffDict_,
            5
        )
    ),
    ct_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ct",
            this->coeffDict_,
            1.87
        )
    ),
    IDDESDelta_(setDelta())
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kOmegaSSTIDDES<BasicTurbulenceModel>::read()
{
    if (kOmegaSSTDES<BasicTurbulenceModel>::read())
    {
        cdt1_.readIfPresent(this->coeffDict());
        cdt2_.readIfPresent(this->coeffDict());
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
