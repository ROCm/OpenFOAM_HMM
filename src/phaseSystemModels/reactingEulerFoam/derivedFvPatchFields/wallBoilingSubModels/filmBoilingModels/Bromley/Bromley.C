/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd
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

#include "Bromley.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace filmBoilingModels
{
    defineTypeNameAndDebug(Bromley, 0);
    addToRunTimeSelectionTable
    (
        filmBoilingModel,
        Bromley,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::filmBoilingModels::Bromley::Bromley
(
    const dictionary& dict
)
:
    filmBoilingModel(),
    Cn_(dict.lookupOrDefault<scalar>("Cn", 0.62))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallBoilingModels::filmBoilingModels::Bromley::~Bromley()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::filmBoilingModels::Bromley::htcFilmBoil
(
    const phaseModel& liquid,
    const phaseModel& vapor,
    const label patchi,
    const scalarField& Tl,
    const scalarField& Tsatw,
    const scalarField& L
) const
{

    const fvPatchScalarField& Tw =
        liquid.thermo().T().boundaryField()[patchi];
    const uniformDimensionedVectorField& g =
        liquid.mesh().time().lookupObject<uniformDimensionedVectorField>("g");

    const scalarField rhoVapor(vapor.thermo().rho(patchi));
    const scalarField rhoLiq(liquid.thermo().rho(patchi));
    const scalarField kappaVapor(vapor.kappa(patchi));

    tmp<volScalarField> tCp = vapor.thermo().Cp();
    const volScalarField& Cp = tCp();
    const scalarField& CpVapor = Cp.boundaryField()[patchi];

    const scalarField muVapor(vapor.mu(patchi));
    const scalarField dbVapor(vapor.d()().boundaryField()[patchi]);

    return
        Cn_*pow
        (
            pow3(kappaVapor)
           *rhoVapor*(rhoLiq - rhoVapor)*mag(g.value())
           *(L + 0.4*CpVapor*max((Tw-Tsatw), scalar(0)))
           /(dbVapor*muVapor*max((Tw-Tsatw), scalar(1e-4))),
            0.25
        );
}


void Foam::wallBoilingModels::filmBoilingModels::Bromley::write
(
    Ostream& os
) const
{
    filmBoilingModel::write(os);
    os.writeKeyword("Cn") << Cn_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //


