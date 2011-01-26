/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

#include "constantHeatTransfer.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(constantHeatTransfer, 0);

addToRunTimeSelectionTable
(
    heatTransferModel,
    constantHeatTransfer,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

constantHeatTransfer::constantHeatTransfer
(
    const surfaceFilmModel& owner,
    const dictionary& dict
)
:
    heatTransferModel(typeName, owner, dict),
    c0_(readScalar(coeffs_.lookup("c0")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

constantHeatTransfer::~constantHeatTransfer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void constantHeatTransfer::correct()
{
    // do nothing
}


tmp<volScalarField> constantHeatTransfer::h() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "htc",
                owner_.time().timeName(),
                owner_.regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            owner_.regionMesh(),
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
