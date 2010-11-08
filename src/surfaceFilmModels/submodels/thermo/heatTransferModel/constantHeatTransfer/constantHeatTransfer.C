/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

#include "constantHeatTransfer.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace surfaceFilmModels
    {
        defineTypeNameAndDebug(constantHeatTransfer, 0);
        addToRunTimeSelectionTable
        (
            heatTransferModel,
            constantHeatTransfer,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceFilmModels::constantHeatTransfer::constantHeatTransfer
(
    const surfaceFilmModel& owner,
    const dictionary& dict
)
:
    heatTransferModel(typeName, owner, dict),
    c0_(readScalar(coeffs_.lookup("c0")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceFilmModels::constantHeatTransfer::~constantHeatTransfer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::surfaceFilmModels::constantHeatTransfer::correct()
{
    // do nothing
}


Foam::tmp<Foam::volScalarField>
Foam::surfaceFilmModels::constantHeatTransfer::h() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "htc",
                owner_.time().timeName(),
                owner_.film(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            owner_.film(),
            dimensionedScalar
            (
                "c0",
                dimEnergy/dimTime/sqr(dimLength)/dimTemperature,
                c0_
            ),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}



// ************************************************************************* //
