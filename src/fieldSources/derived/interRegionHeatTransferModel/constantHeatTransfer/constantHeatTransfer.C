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

#include "constantHeatTransfer.H"
#include "fvm.H"
#include "IObasicSourceList.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcVolumeIntegrate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantHeatTransfer, 0);
    addToRunTimeSelectionTable
    (
        basicSource,
        constantHeatTransfer,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantHeatTransfer::constantHeatTransfer
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    interRegionHeatTransferModel(name, modelType, dict, mesh),
    htCoeffs_(),
    area_()
{
    if (master_)
    {
        htCoeffs_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "htCoeffs",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );

        area_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "area",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );

        htc_.internalField() = htCoeffs_()*area_()/mesh_.V();
        htc_.correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantHeatTransfer::~constantHeatTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::tmp<Foam::volScalarField>Foam::constantHeatTransfer::
calculateHtc()
{
    return htc_;
}


void Foam::constantHeatTransfer::writeData(Ostream& os) const
{
    os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
    interRegionHeatTransferModel::writeData(os);

    os << indent << "constantHeatTransfer";

    dict_.write(os);

    os << decrIndent << indent << token::END_BLOCK << endl;
}


bool Foam::constantHeatTransfer::read(const dictionary& dict)
{
    if (basicSource::read(dict))
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //