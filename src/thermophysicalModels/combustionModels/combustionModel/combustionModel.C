/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

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

#include "combustionModel.H"
#include "surfaceFields.H"
#include "fvScalarMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(combustionModel, 0);
    defineRunTimeSelectionTable(combustionModel, dictionary);
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModel::combustionModel
(
    const dictionary& combustionProps,
    hsCombustionThermo& thermo,
    const compressible::turbulenceModel& turbulence,
    const surfaceScalarField& phi,
    const volScalarField& rho
)
:
    coeffs_(dictionary::null),
    thermo_(thermo),
    turbulence_(turbulence),
    mesh_(phi.mesh()),
    phi_(phi),
    rho_(rho)
{}


Foam::combustionModel::combustionModel
(
    const word& modelType,
    const dictionary& combustionProps,
    hsCombustionThermo& thermo,
    const compressible::turbulenceModel& turbulence,
    const surfaceScalarField& phi,
    const volScalarField& rho
)
:
    coeffs_(combustionProps.subDict(modelType + "Coeffs")),
    thermo_(thermo),
    turbulence_(turbulence),
    mesh_(phi.mesh()),
    phi_(phi),
    rho_(rho)
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::combustionModel::~combustionModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::combustionModel::correct()
{
    // do nothing
}


Foam::tmp<Foam::fvScalarMatrix> Foam::combustionModel::R
(
    volScalarField& Y
) const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix(Y, dimMass/dimTime*Y.dimensions())
    );
}


Foam::tmp<Foam::volScalarField> Foam::combustionModel::dQ() const
{
    return tmp<Foam::volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "dQ",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::combustionModel::wFuelNorm() const
{
    return tmp<Foam::volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "wFuelNorm",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimMass/dimTime/pow3(dimLength), 0.0)
        )
    );
}


bool Foam::combustionModel::read(const dictionary& combustionProps)
{
    coeffs_ = combustionProps.subDict(type() + "Coeffs");

    return true;
}


// ************************************************************************* //
