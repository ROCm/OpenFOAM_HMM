/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "heatTransferCoeff.H"
#include "dictionary.H"
#include "heatTransferCoeffModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(heatTransferCoeff, 0);
    addToRunTimeSelectionTable(functionObject, heatTransferCoeff, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::heatTransferCoeff::calc()
{
    auto& htc = mesh_.lookupObjectRef<volScalarField>(resultName_);

    htcModelPtr_->calc(htc, htcModelPtr_->q());

    htc *= L_/kappa_;

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::heatTransferCoeff::heatTransferCoeff
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict),
    L_(1),
    kappa_(1),
    htcModelPtr_(nullptr)
{
    read(dict);

    setResultName(typeName, "htc:" + htcModelPtr_->type());

    auto* heatTransferCoeffPtr =
        new volScalarField
        (
            IOobject
            (
                resultName_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimPower/dimArea/dimTemperature, Zero)
        );

    mesh_.objectRegistry::store(heatTransferCoeffPtr);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::heatTransferCoeff::read(const dictionary& dict)
{
    if (!fieldExpression::read(dict))
    {
        return false;
    }

    L_ = dict.getCheckOrDefault<scalar>("L", 1, scalarMinMax::ge(0));
    kappa_ =
        dict.getCheckOrDefault<scalar>("kappa", 1, scalarMinMax::ge(SMALL));
    htcModelPtr_ = heatTransferCoeffModel::New(dict, mesh_, fieldName_);

    htcModelPtr_->read(dict);

    return true;
}


// ************************************************************************* //
