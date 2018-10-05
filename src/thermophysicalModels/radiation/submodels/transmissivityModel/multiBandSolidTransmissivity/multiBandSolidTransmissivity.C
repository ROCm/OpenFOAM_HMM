/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenCFD Ltd.
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

#include "multiBandSolidTransmissivity.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(multiBandSolidTransmissivity, 0);

        addToRunTimeSelectionTable
        (
            transmissivityModel,
            multiBandSolidTransmissivity,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::multiBandSolidTransmissivity::multiBandSolidTransmissivity
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    transmissivityModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    tauCoeffs_(),
    nBands_(0)
{
    coeffsDict_.readEntry("transmissivity", tauCoeffs_);
    nBands_ = tauCoeffs_.size();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::multiBandSolidTransmissivity::~multiBandSolidTransmissivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::multiBandSolidTransmissivity::tauEff(const label bandI) const
{
    tmp<volScalarField> tt
    (
        new volScalarField
        (
            IOobject
            (
                "t",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("t", dimless/dimLength, tauCoeffs_[bandI])
        )
    );

    return tt;
}

// ************************************************************************* //
