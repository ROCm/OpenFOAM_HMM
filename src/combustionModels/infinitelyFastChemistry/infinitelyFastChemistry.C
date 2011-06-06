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

#include "infinitelyFastChemistry.H"
#include "addToRunTimeSelectionTable.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{
    defineTypeNameAndDebug(infinitelyFastChemistry, 0);
    addToRunTimeSelectionTable
    (
        combustionModel,
        infinitelyFastChemistry,
        dictionary
    );
};
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::infinitelyFastChemistry::infinitelyFastChemistry
(
    const dictionary& combustionProps,
    hsCombustionThermo& thermo,
    const compressible::turbulenceModel& turbulence,
    const surfaceScalarField& phi,
    const volScalarField& rho
)
:
    combustionModel(typeName, combustionProps, thermo, turbulence, phi, rho),
    C_(readScalar(coeffs_.lookup("C"))),
    singleMixture_
    (
        dynamic_cast<singleStepReactingMixture<gasThermoPhysics>&>(thermo)
    ),
    wFuelNorm_
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
        dimensionedScalar("zero", dimMass/pow3(dimLength)/dimTime, 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::combustionModels::infinitelyFastChemistry::~infinitelyFastChemistry()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::combustionModels::infinitelyFastChemistry::correct()
{
    singleMixture_.fresCorrect();

    const label fuelI = singleMixture_.fuelIndex();

    const volScalarField& YFuel = thermo_.composition().Y()[fuelI];

    const dimensionedScalar s = singleMixture_.s();

    if (thermo_.composition().contains("O2"))
    {
        const volScalarField& YO2 = thermo_.composition().Y("O2");
        wFuelNorm_ == rho_/(mesh_.time().deltaT()*C_)*min(YFuel, YO2/s.value());
    }
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::infinitelyFastChemistry::R(volScalarField& Y) const
{
    const label specieI = thermo_.composition().species()[Y.name()];

    const label fNorm = singleMixture_.specieProd()[specieI];

    const volScalarField fres(singleMixture_.fres(specieI));

    const volScalarField wSpecie
    (
        wFuelNorm_*singleMixture_.specieStoichCoeffs()[specieI]
      / max(fNorm*(Y - fres), scalar(0.001))
    );

    return -fNorm*wSpecie*fres + fNorm*fvm::Sp(wSpecie, Y);
}


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::infinitelyFastChemistry::dQ() const
{
    const label fuelI = singleMixture_.fuelIndex();
    volScalarField& YFuel = thermo_.composition().Y(fuelI);

    return -singleMixture_.qFuel()*(R(YFuel) & YFuel);
}


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::infinitelyFastChemistry::wFuelNorm() const
{
    return wFuelNorm_;
}


bool Foam::combustionModels::infinitelyFastChemistry::read
(
    const dictionary& combustionProps
)
{
    combustionModel::read(combustionProps);
    coeffs_.lookup("C") >> C_ ;

    return true;
}


// ************************************************************************* //
