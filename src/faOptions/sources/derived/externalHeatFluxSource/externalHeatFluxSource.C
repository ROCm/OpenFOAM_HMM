/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "externalHeatFluxSource.H"
#include "fam.H"
#include "faScalarMatrix.H"
#include "physicoChemicalConstants.H"
#include "zeroGradientFaPatchFields.H"
#include "addToRunTimeSelectionTable.H"

using Foam::constant::physicoChemical::sigma;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fa
{
    defineTypeNameAndDebug(externalHeatFluxSource, 0);
    addToRunTimeSelectionTable(option, externalHeatFluxSource, dictionary);
}
}


const Foam::Enum
<
    Foam::fa::externalHeatFluxSource::operationMode
>
Foam::fa::externalHeatFluxSource::operationModeNames
({
    { operationMode::fixedPower, "power" },
    { operationMode::fixedHeatFlux, "flux" },
    { operationMode::fixedHeatTransferCoeff, "coefficient" },
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fa::externalHeatFluxSource::externalHeatFluxSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvPatch& patch
)
:
    fa::faceSetOption(sourceName, modelType, dict, patch),
    mode_(operationModeNames.get("mode", dict)),
    TName_(dict.getOrDefault<word>("T", "T")),
    Q_(0),
    q_(0),
    h_(0),
    Ta_(nullptr),
    emissivity_(dict.getOrDefault<scalar>("emissivity", 0))
{
    fieldNames_.resize(1, TName_);

    fa::option::resetApplied();

    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fa::externalHeatFluxSource::addSup
(
    const areaScalarField& h,
    const areaScalarField& rho,
    faMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (isActive())
    {
        DebugInfo<< name() << ": applying source to "
            << eqn.psi().name() << endl;

        IOobject io
        (
            "Q",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        );

        auto tQ = new areaScalarField
        (
            io,
            regionMesh(),
            dimensionedScalar("q", dimPower/sqr(dimLength), 0),
            zeroGradientFaPatchScalarField::typeName
        );
        areaScalarField& Q = *tQ;

        switch (mode_)
        {
            case fixedPower:
            {
                Q.primitiveFieldRef() = Q_/regionMesh().S().field();
                eqn += Q;

                break;
            }
            case fixedHeatFlux:
            {
                Q.primitiveFieldRef() = q_;
                eqn += Q;
                break;
            }
            case fixedHeatTransferCoeff:
            {
                const dimensionedScalar Ta
                (
                    "Ta",
                    dimTemperature,
                    Ta_->value(mesh_.time().timeOutputValue())
                );

                areaScalarField hp
                (
                    io,
                    regionMesh(),
                    dimensionedScalar
                    (
                        "h",
                        dimPower/sqr(dimLength)/dimTemperature,
                        h_
                    )
                );

                const areaScalarField hpTa(hp*Ta);

                if (emissivity_ > 0)
                {
                    hp -= emissivity_*sigma.value()*pow3(eqn.psi());
                }

                eqn -= fam::SuSp(hp, eqn.psi()) - hpTa;

            }
        }
    }
}


bool Foam::fa::externalHeatFluxSource::read(const dictionary& dict)
{
    if (fa::option::read(dict))
    {
        dict.readIfPresent("T", TName_);
        dict.readIfPresent("emissivity", emissivity_);

        mode_ = operationModeNames.get("mode", dict);

        switch (mode_)
        {
            case fixedPower:
            {
                dict.readEntry("Q", Q_);
                break;
            }
            case fixedHeatFlux:
            {
                dict.readEntry("q", q_);
                break;
            }
            case fixedHeatTransferCoeff:
            {
                dict.readEntry("h", h_);
                Ta_ = Function1<scalar>::New("Ta", dict, &mesh_);
                break;
            }
        }

        return true;
    }

     return false;
}


// ************************************************************************* //
