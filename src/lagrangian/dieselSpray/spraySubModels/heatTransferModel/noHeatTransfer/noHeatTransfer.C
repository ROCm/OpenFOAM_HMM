/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "noHeatTransfer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(noHeatTransfer, 0);

    addToRunTimeSelectionTable
    (
        heatTransferModel,
        noHeatTransfer,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::noHeatTransfer::noHeatTransfer(const dictionary& dict)
:
    heatTransferModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::noHeatTransfer::~noHeatTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::noHeatTransfer::heatTransfer() const
{
    return false;
}


Foam::scalar Foam::noHeatTransfer::Nu(const scalar, const scalar) const
{
    return 0.0;
}


Foam::scalar Foam::noHeatTransfer::relaxationTime
(
    const scalar,
    const scalar,
    const scalar,
    const scalar,
    const scalar,
    const scalar
) const
{
    return GREAT;
}


Foam::scalar Foam::noHeatTransfer::fCorrection(const scalar) const
{
    return 1.0;
}


// ************************************************************************* //
