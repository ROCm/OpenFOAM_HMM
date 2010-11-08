/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

#include "mappedConvectiveHeatTransfer.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "kinematicSingleLayer.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace surfaceFilmModels
    {
        defineTypeNameAndDebug(mappedConvectiveHeatTransfer, 0);
        addToRunTimeSelectionTable
        (
            heatTransferModel,
            mappedConvectiveHeatTransfer,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceFilmModels::mappedConvectiveHeatTransfer::
mappedConvectiveHeatTransfer
(
    const surfaceFilmModel& owner,
    const dictionary& dict
)
:
//    heatTransferModel(typeName, owner, dict),
    heatTransferModel(owner),
    htcConvPrimary_
    (
        IOobject
        (
            "htcConv",
            owner.time().timeName(),
            owner.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        owner.mesh()
    ),
    htcConvFilm_
    (
        IOobject
        (
            htcConvPrimary_.name(), // must have same name as above for mapping
            owner.time().timeName(),
            owner.film(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        owner.film(),
        dimensionedScalar("zero", dimMass/pow3(dimTime)/dimTemperature, 0.0),
        refCast<const kinematicSingleLayer>(owner).pSp().boundaryField().types()
    )
{
    // Update the primary-side convective heat transfer coefficient
    htcConvPrimary_.correctBoundaryConditions();

    // Pull the data from the primary region via direct mapped BCs
    htcConvFilm_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceFilmModels::mappedConvectiveHeatTransfer::
~mappedConvectiveHeatTransfer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::surfaceFilmModels::mappedConvectiveHeatTransfer::correct()
{
    // Update the primary-side convective heat transfer coefficient
    htcConvPrimary_.correctBoundaryConditions();

    // Pull the data from the primary region via direct mapped BCs
    htcConvFilm_.correctBoundaryConditions();
}


Foam::tmp<Foam::volScalarField>
Foam::surfaceFilmModels::mappedConvectiveHeatTransfer::h() const
{
    return htcConvFilm_;
}



// ************************************************************************* //
