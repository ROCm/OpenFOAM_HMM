/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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
    volScalarField& htc = mesh_.lookupObjectRef<volScalarField>(resultName_);

    htcModelPtr_->calc(htc);

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
    htcModelPtr_()
{
    read(dict);

    setResultName(typeName, name + ":htc:" + htcModelPtr_->type());

    volScalarField* heatTransferCoeffPtr
    (
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
            dimensionedScalar("0", dimPower/dimArea/dimTemperature, 0.0)
        )
    );

    mesh_.objectRegistry::store(heatTransferCoeffPtr);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::heatTransferCoeff::~heatTransferCoeff()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::heatTransferCoeff::read(const dictionary& dict)
{
    if (fieldExpression::read(dict))
    {
        htcModelPtr_ = heatTransferCoeffModel::New(dict, mesh_, fieldName_);

        htcModelPtr_->read(dict);

        return true;
    }

    return false;
}


// ************************************************************************* //
