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

#include "tabulatedHeatTransfer.H"

#include "turbulenceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(tabulatedHeatTransfer, 0);
    addToRunTimeSelectionTable
    (
        basicSource,
        tabulatedHeatTransfer,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tabulatedHeatTransfer::tabulatedHeatTransfer
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    interRegionHeatTransferModel(name, modelType, dict, mesh),
    hTable_(),
    area_()
{
    if (master_)
    {
        hTable_.reset
        (
            new interpolation2DTable<scalar>
            (
                dict.subDict(typeName + "Coeffs")
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
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::tabulatedHeatTransfer::~tabulatedHeatTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const Foam::tmp<Foam::volScalarField>Foam::tabulatedHeatTransfer::
calculateHtc()
{
    const fvMesh& secondaryMesh =
        mesh_.time().lookupObject<fvMesh>(mapRegionName());

    const volVectorField& Usecondary =
        secondaryMesh.lookupObject<volVectorField>("U");

    scalarField UMagMapped(htc_.internalField().size(), 0.0);

    secondaryToPrimaryInterpPtr_->interpolateInternalField
    (
        UMagMapped,
        mag(Usecondary),
        meshToMesh::MAP,
        eqOp<scalar>()
    );

    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

    forAll (htc_.internalField(), i)
    {
        htc_.internalField()[i] =
            hTable_->operator()(mag(U[i]), UMagMapped[i]);
    }

    htc_.internalField() = htc_*area_/mesh_.V();

    return htc_;
}


void Foam::tabulatedHeatTransfer::writeData(Ostream& os) const
{
    os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
    interRegionHeatTransferModel::writeData(os);

    os << indent << "tabulatedHeatTransfer";

    dict_.write(os);

    os << decrIndent << indent << token::END_BLOCK << endl;
}


bool Foam::tabulatedHeatTransfer::read(const dictionary& dict)
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