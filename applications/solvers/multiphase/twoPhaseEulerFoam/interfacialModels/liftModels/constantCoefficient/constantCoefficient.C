/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

#include "constantCoefficient.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace liftModels
{
    defineTypeNameAndDebug(constantCoefficient, 0);

    addToRunTimeSelectionTable
    (
        liftModel,
        constantCoefficient,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liftModels::constantCoefficient::constantCoefficient
(
    const dictionary& dict,
    const volScalarField& alpha1,
    const phaseModel& phase1,
    const phaseModel& phase2
)
:
    liftModel(dict, alpha1, phase1, phase2),
    Cl_("Cl", dimless, dict.lookup("Cl"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::liftModels::constantCoefficient::~constantCoefficient()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::liftModels::constantCoefficient::F
(
    const volVectorField& U
) const
{
    return
        Cl_
       *(phase1_*phase1_.rho() + phase2_*phase2_.rho())
       *(
            (phase1_.U() - phase2_.U())
          ^ fvc::curl(U)
        );
}


// ************************************************************************* //
