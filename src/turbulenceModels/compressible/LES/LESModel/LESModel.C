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

#include "LESModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LESModel, 0);
defineRunTimeSelectionTable(LESModel, dictionary);
addToRunTimeSelectionTable(turbulenceModel, LESModel, turbulenceModel);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void LESModel::printCoeffs()
{
    if (printCoeffs_)
    {
        Info<< type() << "Coeffs" << coeffDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * //

LESModel::LESModel
(
    const word& type,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const basicThermo& thermoPhysicalModel
)
:
    turbulenceModel(rho, U, phi, thermoPhysicalModel),

    IOdictionary
    (
        IOobject
        (
            "LESProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    printCoeffs_(lookupOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(subDict(type + "Coeffs")),

    k0_("k0", dimVelocity*dimVelocity, SMALL),

    delta_(LESdelta::New("delta", U.mesh(), *this)),

    wallFunctionDict_(subDict("wallFunctionCoeffs")),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            wallFunctionDict_,
            0.4187
        )
    ),
    E_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "E",
            wallFunctionDict_,
            9.0
        )
    ),
    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            wallFunctionDict_,
            0.07
        )
    ),
    Prt_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Prt",
            wallFunctionDict_,
            0.85
        )
    )
{
    readIfPresent("k0", k0_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void LESModel::correct(const tmp<volTensorField>&)
{
    delta_().correct();
}


void LESModel::correct()
{
    correct(fvc::grad(U_));
}


bool LESModel::read()
{
    if (regIOobject::read())
    {
        coeffDict_ = subDict(type() + "Coeffs");

        readIfPresent("k0", k0_);

        delta_().read(*this);

        wallFunctionDict_ = subDict("wallFunctionCoeffs");
        kappa_.readIfPresent(wallFunctionDict_);
        E_.readIfPresent(wallFunctionDict_);
        Cmu_.readIfPresent(wallFunctionDict_);
        Prt_.readIfPresent(wallFunctionDict_);

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
